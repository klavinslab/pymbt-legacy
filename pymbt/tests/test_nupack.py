from nose.tools import assert_equals
from pymbt import nupack


def test_nupack():
    np_instance = nupack.Nupack(['ATGCGCATGGGAAATAGC', 'ATGCATGCATGCATGC'],
                                material='dna')
    complexes = np_instance.complexes(2)
    concentrations = np_instance.concentrations(2)
    mfe_0 = np_instance.mfe(0)
    mfe_1 = np_instance.mfe(1)
    pairs_0 = np_instance.pairs(0)
    pairs_1 = np_instance.pairs(1)
    np_instance.close()

    assert_equals(complexes['complex_energy'],
                  [-0.49648966, -3.01637728, -10.7753108, -8.59657975,
                   -17.7430411])
    assert_equals(complexes['complexes'],
                  [[1, 0], [0, 1], [2, 0], [1, 1], [0, 2]])
    assert_equals(concentrations['concentration'],
                  [4.671341e-07, 2.733961e-07, 1.642973e-08, 6.386507e-12,
                   1.132988e-07])
    assert_equals(concentrations['energy'],
                  ['-4.964897e-01', '-3.016377e+00', '-1.077531e+01',
                   '-8.596580e+00', '-1.774304e+01'])
    assert_equals(concentrations['types'],
                  [[1, 0], [0, 1], [2, 0], [1, 1], [0, 2]])
    assert_equals(mfe_0, 0.0)
    assert_equals(mfe_1, -2.379)
    assert_equals(pairs_0['probabilities'],
                  [0.95587, 0.90174, 0.92009, 0.79434, 0.9299, 0.86772,
                   0.86381, 0.79284, 0.96115, 0.91567, 0.91735, 0.95318,
                   0.93185, 0.89207, 0.87689, 0.99119, 0.84085, 0.83946])
    assert_equals(pairs_0['type'],
                  [(1, 19), (2, 19), (3, 19), (4, 19), (5, 19), (6, 19),
                   (7, 19), (8, 19), (9, 19), (10, 19), (11, 19), (12, 19),
                   (13, 19), (14, 19), (15, 19), (16, 19), (17, 19), (18, 19)])
    assert_equals(pairs_1['probabilities'],
                  [0.6304, 0.55413, 0.080166, 0.05002, 0.092026, 0.5513,
                   0.55794, 0.97272, 0.97456, 0.54112, 0.50081, 0.080504,
                   0.11048, 0.16068, 0.52622, 0.56212])
    assert_equals(pairs_1['type'],
                  [(1, 17), (2, 17), (3, 17), (4, 17), (5, 17), (6, 17),
                  (7, 17), (8, 17), (9, 17), (10, 17), (11, 17), (12, 17),
                  (13, 17), (14, 17), (15, 17), (16, 17)])
