"""RTN basis: orthonormality, handedness, known geometry, covariance rotation (spec sections 12, 36)."""

import numpy as np
import pytest

from tart_gnss.stage6.errors import RTNGeometryError
from tart_gnss.stage6.rtn import covariance_to_rtn, rtn_basis, to_rtn


def test_orthonormal_right_handed():
    rng = np.random.default_rng(0)
    for _ in range(50):
        r, v = rng.normal(size=3) * 2e7, rng.normal(size=3) * 3e3
        q = rtn_basis(r, v)
        assert np.allclose(q @ q.T, np.eye(3), atol=1e-12)
        assert np.linalg.det(q) == pytest.approx(1.0, abs=1e-12)
        assert np.allclose(np.cross(q[0], q[1]), q[2], atol=1e-12)
        assert to_rtn(q, r)[0] == pytest.approx(np.linalg.norm(r))
        assert abs(to_rtn(q, v)[2]) < 1e-6 and to_rtn(q, v)[1] > 0  # velocity has no N part, positive T


def test_known_circular_orbit():
    q = rtn_basis(np.array([7e6, 0, 0]), np.array([0, 7.5e3, 0]))
    assert np.allclose(q, np.eye(3))


def test_degenerate_raises():
    with pytest.raises(RTNGeometryError):
        rtn_basis(np.array([7e6, 0, 0]), np.array([1.0, 0, 0]))
    with pytest.raises(RTNGeometryError):
        rtn_basis(np.zeros(3), np.array([0, 1.0, 0]))


def test_covariance_rotation_numerically():  # 34.5
    rng = np.random.default_rng(1)
    q = rtn_basis(rng.normal(size=3) * 2e7, rng.normal(size=3) * 3e3)
    a = rng.normal(size=(3, 3))
    c = a @ a.T * 1e6
    x = rng.multivariate_normal(np.zeros(3), c, size=40000)
    emp = np.cov((x @ q.T).T)
    rot = covariance_to_rtn(q, c)
    assert np.allclose(emp, rot, rtol=0.05, atol=0.03 * np.max(np.abs(rot)))
    assert np.allclose(np.linalg.eigvalsh(rot), np.linalg.eigvalsh(c), rtol=1e-10)
