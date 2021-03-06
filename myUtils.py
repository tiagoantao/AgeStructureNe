try:
    import configparser as ConfigParser
except:
    import ConfigParser
from functools import reduce

meta = {
    "bighorn": {'name': 'Bighorn sheep'},
    "bison": {'name': ''},
    "bottlenose": {'name': 'Bottlenose dolphin'},
    "btrout": {'name': ''},
    "bullt2": {'name': ''},
    "bulltrout": {'name': ''},
    "bullpred": {'name': ''},
    "cod": {'name': ''},
    "emperor": {'name': ''},
    "fraley": {'name': ''},
    "grizzly": {'name': 'Grizzly bear'},
    "lake": {'name': 'Lake trout'},
    "lturtle": {'name': ''},
    "mcrab": {'name': ''},
    "mosquito": {'name': ''},
    "primrose": {'name': ''},
    "sagebrush": {'name': ''},
    "sagegrouse": {'name': ''},
    "sardine": {'name': 'Sardine'},
    "seaweed": {'name': ''},
    "shepard": {'name': ''},
    "snake": {'name': ''},
    "sparrow": {'name': ''},
    "surchin": {'name': 'Sea urchin'},
    "wfrog": {'name': ''},
}

ages = {
    "bighorn": 16,
    "bison": 20,
    "bottlenose": 36,
    "btrout": 15,
    "bullt2": 22,
    "bulltrout": 11,
    "bullpred": 11,
    "cod": 21,
    "emperor": 26,
    "fraley": 11,
    "grizzly": 32,
    "lake": 63,
    "lturtle": 55,
    "mcrab": 20,
    "mosquito": 32,
    "primrose": 23,
    "sagebrush": 19,
    "sagegrouse": 16,
    "sardine": 14,
    "seaweed": 23,
    "synseaweed": 23,
    "shepard": 11,
    "snake": 9,
    "sparrow": 7,
    "surchin": 6,
    "wfrog": 4,
}

gammaAFemale = {
    "bison": {
        1: [None] * 19
    },
    "ecology": {
        3: [0.318, 0.636, 0.954, 1.271]
    },
    "sardine": {
        4: [9.000, 0.113, 0.299, 0.462, 0.655, 0.900, 1.147, 1.248,
            1.456, 1.663, 1.663, 1.665, 1.661]
    },
    "sparrow": {
        1: [None] * 6,
        2: [9.0, 5.090, 5.505, 5.839, 6.260, 6.680]
    }
}

gammaBFemale = {
    "bison": {
        1: [None] * 19
    },
    "ecology": {
        3: [2.0] * 4
    },
    "sardine": {
        4: [3.0] * 13,
    },
    "sparrow": {
        1: [None] * 6,
        2: [1.0] * 6,
    },
}

gammaAMale = {
    "bison": {
        4: [9.0, 9.0, 9.0, 0.046, 0.152, 0.322, 0.461, 0.562, 0.632,
            0.670, 0.676, 0.639, 0.575, 0.480, 0.349, 0.179, 9.000,
            9.000, 9.000]
    },
    "ecology": {
        3: [0.136, 0.272, 0.409, 0.545]
    },
    "sardine": {
        8: [9.000, 0.048, 0.128, 0.198, 0.281, 0.386, 0.492, 0.535,
            0.624, 0.713, 0.713, 0.714, 0.712]

    },
    "sparrow": {
        2: [9.000, 5.090, 5.505, 5.839, 6.260, 6.680],
        4: [9.000, 1.697, 1.835, 1.946, 2.087, 2.227],
    }
}

gammaBMale = {
    "bison": {
        4: [3.0] * 19
    },
    "ecology": {
        3: [2.0] * 4
    },
    "sardine": {
        8: [7.0] * 13,
    },
    "sparrow": {
        2: [1.0] * 6,
        4: [3.0] * 6,
    },
}


survivalFemale = {
    "bighorn": [0.765, 0.792, 0.792, 0.792, 0.792, 0.792, 0.792, 0.735,
                0.735, 0.735, 0.735, 0.735, 0.484, 0.484, 0.000],
    "bison": [0.804, 0.819, 0.829, 0.835, 0.837, 0.836, 0.831, 0.822, 0.809,
              0.791, 0.767, 0.736, 0.698, 0.65, 0.593, 0.527, 0, 0],
    "bottlenose": [0.842, 0.88, 0.909, 0.93, 0.945, 0.955, 0.962, 0.966,
                   0.968, 0.969, 0.968, 0.967, 0.964, 0.96, 0.956, 0.951,
                   0.945, 0.937, 0.929, 0.92, 0.909, 0.898, 0.884, 0.869,
                   0.852, 0.833, 0.812, 0.789, 0.763, 0.735, 0.704, 0.67,
                   0.63, 0.594],
    "btrout": [0.48]*5+[0.5]*8,
    "bullt2": [0.44, 0.21, 0.21, 0.15, 0.45, 0.5, 0.55, 0.45, 0.53, 0.6,
               0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.85, 0.7, 0.6, 0.2],
    "bulltrout": [0.48, 0.45, 0.447, 0.453, 0.462, 0.470, 0.477, 0.482, 0.485],
    "bullpred": [0.48, 0.45, 0.447, 0.225, 0.225, 0.470, 0.477, 0.482, 0.485],
    "cod": [0.1, 0.5] + [0.82] * 17,
    "emperor": [0.86] * 24,
    "fraley": [0.22, 0.32, 0.39, 0.399, 0.388, 0.377, 0.368, 0.362, 0.358],
    "grizzly": [0.817] + 29*[0.922],
    "lake": [0.45, 0.78] + 59*[0.92],
    "lturtle": [0.786, 0.786, 0.834, 0.825, 0.818, 0.812, 0.808, 0.804,
                0.801, 0.799, 0.760, 0.676, 0.678, 0.679, 0.681, 0.683,
                0.685, 0.687, 0.743, 0.747, 0.752, 0.756, 0.761, 0.765,
                0.770, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809,
                0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809,
                0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809,
                0.809, 0.809, 0.809, 0.809, 0.809],
    "mcrab": 18*[0.74],
    "mosquito": 19*[0.85] + 4*[0.90] + 4*[0.62] + 3*[0.37],
    "primrose": [0.492, 0.605, 0.683, 0.734, 0.765, 0.784, 0.794, 0.801,
                 0.804, 0.806, 0.808, 0.808, 0.809, 0.809, 0.809, 0.809, 0.810,
                 0.810, 0.810, 0.810, 0.810],
    "sagebrush": [0.507, 0.507, 0.778, 0.926, 0.926] + [0.941]*12,
    "sagegrouse": [0.731] + [0.733]*13,
    "sardine": [0.47]+[0.67]*11,
    "seaweed": [0.63, 0.674, 0.708, 0.734, 0.756, 0.773, 0.787, 0.798, 0.805,
                0.811, 0.815, 0.818, 0.82, 0.822, 0.823, 0.823, 0.824, 0.824,
                0.825, 0.825, 0.825],
    "synseaweed": [0.1, 0.674, 0.708, 0.734, 0.756, 0.773, 0.787, 0.798,
                   0.805, 0.811, 0.815, 0.818, 0.82, 0.822, 0.823, 0.823,
                   0.824, 0.824, 0.825, 0.825, 0.825],
    "shepard": [0.22, 0.32, 0.39, 0.399, 0.388, 0.377, 0.368, 0.362, 0.358],
    "snake": [0.486, 0.485, 0.480, 0.500, 0.500, 0.500, 0.333],
    "sparrow": [0.18, 0.528, 0.537, 0.529, 0.519],
    "surchin": [0.323, 0.231, 0.176, 0.141],
    "wfrog": [0.143, 0.124],
}

survivalMale = {
    "bighorn": [0.826, 0.850, 0.850, 0.850, 0.850, 0.850, 0.850, 0.598,
                0.598, 0.598, 0.598, 0.598, 0.000, 0.000, 0.000],
    "bison": [0.796, 0.791, 0.785, 0.78, 0.774, 0.768, 0.762, 0.756, 0.75,
              0.744, 0.737, 0.731, 0.725, 0.718, 0.711, 0.704, 0.697, 0.69],
    "bottlenose": [0.824, 0.851, 0.873, 0.89, 0.903, 0.913, 0.919, 0.924,
                   0.926, 0.927, 0.926, 0.924, 0.921, 0.916, 0.911, 0.904,
                   0.896, 0.887, 0.877, 0.867, 0.854, 0.841, 0.827, 0.811,
                   0.794, 0.776, 0.756, 0.735, 0.712, 0.687, 0.661, 0.634,
                   0, 0],
    "btrout": [0.48]*5+[0.5]*8,
    "bullt2": [0.44, 0.21, 0.21, 0.15, 0.45, 0.5, 0.55, 0.45, 0.53, 0.6,
               0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.85, 0.7, 0.6, 0.2],
    "bulltrout": [0.48, 0.45, 0.447, 0.453, 0.462, 0.470, 0.477, 0.482, 0.485],
    "bullpred": [0.48, 0.45, 0.447, 0.225, 0.225, 0.470, 0.477, 0.482, 0.485],
    "cod": [0.1, 0.5] + [0.82] * 17,
    "emperor": [0.86] * 24,
    "fraley": [0.22, 0.32, 0.39, 0.399, 0.388, 0.377, 0.368, 0.362, 0.358],
    "grizzly": [0.817] + 25 * [0.823] + 4 * [0],
    "lake": [0.45, 0.78] + 59*[0.92],
    "lturtle": [0.786, 0.786, 0.834, 0.825, 0.818, 0.812, 0.808, 0.804,
                0.801, 0.799, 0.760, 0.676, 0.678, 0.679, 0.681, 0.683,
                0.685, 0.687, 0.743, 0.747, 0.752, 0.756, 0.761, 0.765,
                0.770, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809,
                0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809,
                0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809, 0.809,
                0.809, 0.809, 0.809, 0.809, 0.809],
    "mcrab": 18*[0.68],
    "mosquito": 19*[0.85] + 4*[0.90] + 4*[0.62] + 3*[0.37],
    "primrose": [0.492, 0.605, 0.683, 0.734, 0.765, 0.784, 0.794, 0.801,
                 0.804, 0.806, 0.808, 0.808, 0.809, 0.809, 0.809, 0.809, 0.810,
                 0.810, 0.810, 0.810, 0.810],
    "sagebrush": [0.507, 0.507, 0.778, 0.926, 0.926] + [0.941]*12,
    "sagegrouse": [0.731] + [0.733]*13,
    "sardine": [0.47]+[0.67]*11,
    "seaweed": [0.63, 0.674, 0.708, 0.734, 0.756, 0.773, 0.787, 0.798,
                0.805, 0.811, 0.815, 0.818, 0.82, 0.822, 0.823, 0.823, 0.824,
                0.824, 0.825, 0.825, 0.825],
    "synseaweed": [0.10, 0.674, 0.708, 0.734, 0.756, 0.773, 0.787, 0.798,
                   0.805, 0.811, 0.815, 0.818, 0.82, 0.822, 0.823, 0.823,
                   0.824, 0.824, 0.825, 0.825, 0.825],
    "shepard": [0.22, 0.32, 0.39, 0.399, 0.388, 0.377, 0.368, 0.362, 0.358],
    "snake": [0.486, 0.485, 0.480, 0.500, 0.500, 0.500, 0.333],
    "sparrow": [0.18, 0.528, 0.537, 0.529, 0.519],
    "surchin": [0.323, 0.231, 0.176, 0.141],
    "wfrog": [0.143, 0.124],
}

fecFemale = {
    "bighorn": [0.000, 0.073, 0.218, 0.334, 0.334, 0.334, 0.334, 0.334, 0.334,
                0.334, 0.334, 0.334, 0.331, 0.331, 0.331],
    "bison": [0, 0.08, 0.75, 0.93, 0.93, 0.93, 0.93, 0.93, 0.92, 0.92, 0.92,
              0.92, 0.92, 0.8, 0.8, 0.5, 0.5, 0, 0],
    "bottlenose": [0, 0, 0, 0, 0, 0, 0, 0.172, 0.172, 0.171, 0.171, 0.171,
                   0.142, 0.132, 0.135, 0.134, 0.135, 0.133, 0.133, 0.131,
                   0.132, 0.13, 0.129, 0.129, 0.128, 0.127, 0.126, 0.124,
                   0.124, 0.123, 0.121, 0.121, 0.118, 0.117, 0.116],
    "btrout": [0]*5 + [1.9, 26.8, 43.5, 65.5, 73.7] + [72.7]*4,
    "bullt2": [0.0, 0.0, 0.0, 0.156, 0.352, 0.536, 0.686, 0.796, 0.872, 0.921,
               0.953, 0.96, 0.97, 0.98, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,
               0.99],
    "bulltrout": [0.0, 0.0, 0.156, 0.352, 0.536, 0.686, 0.796, 0.872, 0.921,
                  0.953],
    "bullpred": [0.0, 0.0, 0.156, 0.352, 0.536, 0.686, 0.796, 0.872, 0.921,
                 0.953],
    "cod": [0, 0, 0, 0, 0, 0, 1.049, 1.275, 1.543, 1.817, 1.974, 2.425,
            2.883, 3.438, 4.237, 4.676, 6.007, 6.486, 6.966, 7.499],
    "emperor": [0, 0, 0.034, 0.121, 0.205, 0.223, 0.249, 0.272, 0.293, 0.311,
                0.326, 0.340, 0.353, 0.363, 0.373, 0.381, 0.389, 0.395, 0.401,
                0.406, 0.410, 0.414, 0.417, 0.420, 0.423],
    "fraley": [0, 0, 265, 589, 964, 1344, 1555, 1771, 1963, 1966],
    "grizzly": 3*[0] + 28*[0.318],
    "lake": 6*[0] + [1058, 5381, 10367, 9789] + 52*[12793],
    "lturtle": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.9, 8.5, 14.3, 20.0, 25.6,
                30.9, 103.2, 36.7, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0,
                80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0,
                80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0],
    "mcrab": [0, 0, 0, 0, 0, 0, 0, 0, 2.772, 3.463, 4.012, 6.918, 8.694,
              10.926, 13.734, 15.996, 15.996, 15.996, 18.628],
    "mosquito": 19*[0] + 4*[0.2] + 4*[0.171] + 4*[0.146],
    "primrose": [0.000, 0.004, 0.023, 0.058, 0.095, 0.126, 0.148, 0.162,
                 0.171, 0.177, 0.180, 0.182, 0.184, 0.184, 0.185, 0.185, 0.185,
                 0.185, 0.185, 0.185, 0.186, 0.186],
    "sagebrush": [11, 17, 46, 57, 75, 77, 80, 83, 86, 88, 82, 72, 66, 62, 58,
                  54, 29, 0],
    "sagegrouse": [0.55] + [0.817]*14,
    "sardine": [0, 146, 388, 599, 849, 1167, 1487, 1617, 1887, 2156, 2156,
                2156, 2156],
    "seaweed": [0, 0.05, 0.17, 0.34, 0.51, 0.64, 0.75, 0.83, 0.89, 0.94, 0.97,
                0.99, 1.01, 1.02, 1.03, 1.04, 1.04, 1.04, 1.04, 1.05, 1.05,
                1.05],
    "synseaweed": [0, 0.05, 0.17, 0.34, 0.51, 0.64, 0.75, 0.83, 0.89, 0.94,
                   0.97, 0.99, 1.01, 1.02, 1.03, 1.04, 1.04, 1.04, 1.04, 1.05,
                   1.05, 1.05],
    "shepard": [0, 0, 7.5, 14.9, 19.2, 21.6, 22.9, 23.7, 24.2, 24.5],
    "snake": [1.03, 2.81, 3.25, 3.52, 3.67, 3.75, 3.81, 3.84],
    "sparrow": [0, 2.546, 2.754, 2.921, 3.13, 3.339],
    "surchin": [2.64, 4.47, 4.88, 4.97, 4.99],
    "wfrog": [1.04, 286.0, 342.0],
}

fecMale = {
    "bighorn": [0.000, 0.073, 0.218, 0.334, 0.334, 0.334, 0.334, 0.334, 0.334,
                0.334, 0.334, 0.334, 0.331, 0.000, 0.000],
    "bison": [0, 0, 0, 0.11, 0.36, 0.76, 1.09, 1.33, 1.5, 1.59, 1.59, 1.52,
              1.36, 1.13, 0.82, 0.42, 0, 0, 0],
    "bottlenose": [0, 0, 0, 0, 0, 0, 0, 0.172, 0.172, 0.171, 0.171, 0.171,
                   0.142, 0.132, 0.135, 0.134, 0.135, 0.133, 0.133, 0.131,
                   0.132, 0.13, 0.129, 0.129, 0.128, 0.127, 0.126, 0.124,
                   0.124, 0.123, 0.121, 0.121, 0.118, 0.117, 0.116],
    "btrout": [0] * 4 + [1.9] + [11.7, 17.8, 26.7, 34.4, 40.7] + [67.6] * 4,
    "bullt2": [0.0, 0.0, 0.0, 0.156, 0.352, 0.536, 0.686, 0.796, 0.872, 0.921,
               0.953, 0.96, 0.97, 0.98, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,
               0.99],
    "bulltrout": [0.0, 0.0, 0.156, 0.352, 0.536, 0.686, 0.796, 0.872, 0.921,
                  0.953],
    "bullpred": [0.0, 0.0, 0.156, 0.352, 0.536, 0.686, 0.796, 0.872, 0.921,
                 0.953],
    "cod": [0, 0, 0, 0, 0, 0, 1.049, 1.275, 1.543, 1.817, 1.974, 2.425,
            2.883, 3.438, 4.237, 4.676, 6.007, 6.486, 6.966, 7.499],
    "emperor": [0, 0, 0.034, 0.121, 0.205, 0.223, 0.249, 0.272, 0.293, 0.311,
                0.326, 0.340, 0.353, 0.363, 0.373, 0.381, 0.389, 0.395, 0.401,
                0.406, 0.410, 0.414, 0.417, 0.420, 0.423],
    "fraley": [0, 0, 265, 589, 964, 1344, 1555, 1771, 1963, 1966],
    "grizzly": 3 * [0] + 24 * [0.318] + 4 * [0],
    "lake": 6*[0] + [1058, 5381, 10367, 9789] + 52*[12793],
    "lturtle": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.9, 8.5, 14.3, 20.0, 25.6,
                30.9, 103.2, 36.7, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0,
                80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0,
                80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0, 80.0],
    "mcrab": [0, 0, 0, 0, 0, 0, 0, 0, 2.772, 3.463, 4.012, 6.918, 8.694,
              10.926, 13.734, 15.996, 15.996, 15.996, 18.628],
    "mosquito": 19*[0] + 4*[0.2] + 4*[0.171] + 4*[0.146],
    "primrose": [0.000, 0.004, 0.023, 0.058, 0.095, 0.126, 0.148, 0.162,
                 0.171, 0.177, 0.180, 0.182, 0.184, 0.184, 0.185, 0.185, 0.185,
                 0.185, 0.185, 0.185, 0.186, 0.186],
    "sagebrush": [11, 17, 46, 57, 75, 77, 80, 83, 86, 88, 82, 72, 66, 62, 58,
                  54, 29, 0],
    "sagegrouse": [0.55] + [0.817]*14,
    "sardine": [0, 146, 388, 599, 849, 1167, 1487, 1617, 1887, 2156, 2156,
                2156, 2156],
    "seaweed": [0, 0.05, 0.17, 0.34, 0.51, 0.64, 0.75, 0.83, 0.89, 0.94, 0.97,
                0.99, 1.01, 1.02, 1.03, 1.04, 1.04, 1.04, 1.04, 1.05, 1.05,
                1.05],
    "synseaweed": [0, 0.05, 0.17, 0.34, 0.51, 0.64, 0.75, 0.83, 0.89, 0.94,
                   0.97, 0.99, 1.01, 1.02, 1.03, 1.04, 1.04, 1.04, 1.04, 1.05,
                   1.05, 1.05],
    "shepard": [0, 0, 7.5, 14.9, 19.2, 21.6, 22.9, 23.7, 24.2, 24.5],
    "snake": [1.03, 2.81, 3.25, 3.52, 3.67, 3.75, 3.81, 3.84],
    "sparrow": [0, 2.546, 2.754, 2.921, 3.13, 3.339],
    "surchin": [2.64, 4.47, 4.88, 4.97, 4.99],
    "wfrog": [1.04, 286.0, 342.0],
}


def getVarConf(fName):
    config = ConfigParser.ConfigParser()
    config.read(fName)

    N0 = eval(config.get("var", "N0"))
    sampCohort = eval(config.get("var", "sampCohort"))
    sampSize = eval(config.get("var", "sampSize"))
    try:
        sampSNP = eval(config.get("var", "sampSNP"))
    except:
        sampSNP = {}

    numGens = config.getint("sim", "numGens")
    reps = config.getint("sim", "reps")
    dataDir = config.get("sim", "dataDir")

    return N0, sampCohort, sampSize, sampSNP, numGens, reps, dataDir


def getConfig(fName):
    config = ConfigParser.ConfigParser()
    config.read(fName)

    class Cfg:
        pass

    cfg = Cfg()

    cfg.popSize = config.getint("pop", "popSize")
    cfg.ages = eval(config.get("pop", "ages"))
    if config.has_option("pop", "isMonog"):
        cfg.isMonog = config.getboolean("pop", "isMonog")
    else:
        cfg.isMonog = False
    if config.has_option("pop", "forceSkip"):
        cfg.forceSkip = config.getfloat("pop", "forceSkip") / 100
    else:
        cfg.forceSkip = 0
    if config.has_option("pop", "skip"):
        cfg.skip = eval(config.get("pop", "skip"))
    else:
        cfg.skip = []
    if config.has_option("pop", "litter"):
        cfg.litter = eval(config.get("pop", "litter"))
    else:
        cfg.litter = None
    if config.has_option("pop", "male.probability"):
        cfg.maleProb = config.getfloat("pop", "male.probability")
    else:
        cfg.maleProb = 0.5
    if config.has_option("pop", "gamma.b.male"):
        cfg.doNegBinom = True
        cfg.gammaAMale = eval(config.get("pop", "gamma.a.male"))
        cfg.gammaBMale = eval(config.get("pop", "gamma.b.male"))
        cfg.gammaAFemale = eval(config.get("pop", "gamma.a.female"))
        cfg.gammaBFemale = eval(config.get("pop", "gamma.b.female"))
    else:
        cfg.doNegBinom = False
    cfg.N0 = config.getint("pop", "N0")
    if config.has_option("pop", "survival"):
        cfg.survivalMale = eval(config.get("pop", "survival"))
        cfg.survivalFemale = eval(config.get("pop", "survival"))
    else:
        cfg.survivalMale = eval(config.get("pop", "survival.male"))
        cfg.survivalFemale = eval(config.get("pop", "survival.female"))
    cfg.fecundityMale = eval(config.get("pop", "fecundity.male"))
    cfg.fecundityFemale = eval(config.get("pop", "fecundity.female"))
    if config.has_option("pop", "startLambda"):
        cfg.startLambda = config.getint("pop", "startLambda")
        cfg.lbd = config.getfloat("pop", "lambda")
        # cfg.lbd = mp.mpf(config.get("pop", "lambda"))
    else:
        cfg.startLambda = 99999
        # cfg.lbd = mp.mpf(1.0)
        cfg.lbd = 1.0
    if config.has_option("pop", "Nb"):
        cfg.Nb = config.getint("pop", "Nb")
        cfg.NbVar = config.getfloat("pop", "NbVar")
    else:
        cfg.Nb = None
        cfg.NbVar = None

    cfg.startAlleles = config.getint("genome", "startAlleles")
    cfg.mutFreq = config.getfloat("genome", "mutFreq")
    cfg.numMSats = config.getint("genome", "numMSats")
    if config.has_option("genome", "numSNPs"):
        cfg.numSNPs = config.getint("genome", "numSNPs")
    else:
        cfg.numSNPs = 0

    cfg.reps = config.getint("sim", "reps")
    if config.has_option("sim", "startSave"):
        cfg.startSave = config.getint("sim", "startSave")
    else:
        cfg.startSave = 0
    cfg.gens = config.getint("sim", "gens")
    cfg.dataDir = config.get("sim", "dataDir")

    return cfg

def getAgeFecund(cfg):
    ages = []
    for fec in [cfg.fecundityMale, cfg.fecundityFemale]:
        for i in range(len(fec)):
            if fec[i]>0:
                ages.append(i+1)
                break
    return tuple(ages)


def calcHarmonic(lst):
    #Negatives are converted to infs (LDNe)
    myLst = [x < 0 and float("inf") or x for x in lst]
    try:
        return 1.0 * len(lst) / reduce(lambda x, y: 1.0 / y + x, myLst, 0.0)
    except ZeroDivisionError:
        return None


def getNonGenStats(f):
    l = f.readline()
    stats = []
    while l != '':
        toks = l.rstrip().split(' ')
        try:
            gender = int(toks[3])
            stats.append({
                "gen": int(toks[0]),
                "rep": int(toks[1]),
                "id": int(float(toks[2])),
                "gender": gender,
                "father": int(float(toks[4])),
                "mother": int(float(toks[5])),
                "age": int(float(toks[6])),
            })
        except IndexError:
            break

        l = f.readline()
    return stats


def getDemo(project, model, N0, rep):
    f = open('data/%s/%d%s%d.demo' % (project, N0, model, rep))
    for l in f:
        toks = [x for x in l.strip().rstrip().split(' ') if x != '']
        cycle = int(toks[0])
        gentime = float(toks[1])
        pyramid = [int(x) for x in toks[2:]]
        yield {'cycle': cycle, 'generation': gentime, 'pyramid': pyramid}
    f.close()


def getVk(project, model, N0, rep):
    f = open('data/%s/%d%s%d.vk' % (project, N0, model, rep))
    for l in f:
        toks = [x for x in l.strip().rstrip().split(' ') if x != '']
        cycle = int(toks[0])
        k = float(toks[1])
        vk = float(toks[2])
        nb = float(toks[3])
        ne = float(toks[4])
        yield {'cycle': cycle, 'k': k, 'vk': vk, 'nb': nb, 'ne': ne}
    f.close()
