import numpy as np

from lib.vector import Vector

def test_vector_angle():
    v1 = Vector(1, 0)
    v2 = Vector(0, 1)
    v3 = Vector(-1, 0)
    v4 = Vector(0, -1)

    assert v1.angle == np.deg2rad(0)
    assert v2.angle == np.deg2rad(90)
    assert v3.angle == np.deg2rad(180)
    assert v4.angle == np.deg2rad(-90)
def test_vector_perturbation():
    v1 = Vector(1, 0)
    v2 = Vector(0, 1)
    v3 = Vector(-1, 0)
    v4 = Vector(0, -1)

    v1_perturbed = v1.perturbate((np.deg2rad(-15), np.deg2rad(15)))
    v2_perturbed = v2.perturbate((np.deg2rad(-15), np.deg2rad(15)))
    v3_perturbed = v3.perturbate((np.deg2rad(-15), np.deg2rad(15)))
    v4_perturbed = v4.perturbate((np.deg2rad(-15), np.deg2rad(15)))

    assert np.rad2deg(v1_perturbed.angle) > -15 and np.rad2deg(v1_perturbed.angle) < 15
    assert np.rad2deg(v2_perturbed.angle) > 75 and np.rad2deg(v2_perturbed.angle) < 105
    assert np.rad2deg(v3_perturbed.angle) > 165 and np.rad2deg(v3_perturbed.angle) < 195
    assert np.rad2deg(v4_perturbed.angle) > -105 and np.rad2deg(v4_perturbed.angle) < -75