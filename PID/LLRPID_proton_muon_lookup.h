#ifndef PROTON_MUON_LOOKUP_H
#define PROTON_MUON_LOOKUP_H
#include <stdlib.h>
#include <vector>

namespace searchingfornues
{
  struct ProtonMuonLookUpParameters
  {

    std::vector<float> dedx_edges_pl_0 = {
      0.000, 0.500, 1.000, 1.500, 2.000, 2.500, 3.000, 3.500, 4.000, 4.500, 5.000, 5.500, 6.000, 6.500, 7.000, 7.500, 8.000, 9.000, 10.000, 12.000,
      15.000, 20.000, 25.000, 30.000, 35.000, 40.000, 45.000, 50.000,
    };



    std::vector<std::vector<float>> parameters_edges_pl_0 = {
      {0.000, 2.000, 4.000, 7.000, 10.000, 15.000, 20.000, 30.000, 50.000, 100.000, 300.000, 2000.000,},
      {0.000, 0.600, 1.000, 1.500, 3.000, 30.000,}
    };


    std::vector<float> dedx_pdf_pl_0 = {
      0.806, 1.665, 2.779, 2.624, 2.178, 1.715, 1.461, 1.547, 1.594, 1.644, 1.655, 1.567, 1.426, 1.310, 1.100, 1.012, 0.687, 0.160, -0.652, -1.623,
      -2.193, -2.402, -2.195, -1.919, -1.055, -0.869, -0.030, 1.389, 1.488, 2.439, 2.410, 2.069, 1.662, 1.559, 1.388, 1.454, 1.185, 1.235, 1.067, 0.814,
      0.506, 0.423, 0.036, -0.411, -0.701, -1.447, -1.744, -2.105, -2.107, -1.889, -1.305, -0.886, -0.886, -0.541, 1.468, 1.723, 1.816, 1.635, 1.688, 1.704,
      1.069, 0.676, 0.798, 0.433, 0.300, 0.139, 0.213, -0.405, -0.371, -0.592, -1.096, -1.228, -1.639, -1.853, -2.187, -2.275, -2.397, -3.362, -1.164, -2.715,
      0.000, 0.829, 0.970, 0.839, 1.086, 0.858, 0.443, -0.022, -0.051, -0.125, -0.711, -1.010, -1.533, -1.044, -2.467, -1.572, -1.790, -1.847, -2.324, -2.160,
      -2.285, -2.653, -3.598, 0.000, 0.000, 0.000, 0.000, 0.000, 0.180, -0.019, -0.150, -0.306, -0.682, -1.011, -1.987, -1.757, -2.641, -1.437, 0.000, 0.000,
      -2.255, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -1.617, 1.915, 3.145, 3.144, 2.312,
      2.056, 2.143, 2.182, 2.123, 1.732, 1.436, 0.772, 0.293, -0.291, -0.966, -1.528, -2.448, -3.300, -3.685, -3.455, -2.563, -1.452, -0.856, -0.906, -1.530,
      -1.064, -0.049, 1.688, 1.643, 2.474, 2.754, 2.207, 1.971, 1.897, 1.829, 1.551, 1.110, 0.923, 0.449, 0.123, -0.282, -0.990, -1.353, -2.367, -2.850,
      -3.008, -2.884, -2.621, -1.825, -1.157, -1.682, -1.921, -0.989, -0.439, 1.442, 1.783, 2.320, 2.163, 2.393, 1.805, 1.395, 1.140, 0.990, 0.572, 0.263,
      0.204, -0.305, -0.976, -1.258, -1.728, -2.013, -2.539, -2.935, -2.739, -2.899, -3.362, -2.693, -2.133, -0.054, -0.747, 0.000, 0.821, 1.145, 1.259, 1.372,
      1.279, 1.194, 0.407, 0.387, 0.082, -0.311, -0.614, -1.127, -1.505, -1.302, -2.010, -2.620, -2.450, -2.514, -2.934, -3.490, -3.316, -4.069, -3.392, 0.000,
      0.000, 0.000, 0.000, 0.181, -0.176, -0.197, -0.156, -0.138, -0.883, -1.799, -0.595, -1.288, -1.442, 0.000, -2.441, 0.000, -2.061, 0.000, 0.000, -1.981,
      0.000, -1.799, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.108, 2.494, 3.166, 3.132, 2.762, 2.726, 2.529, 2.067, 1.405, 0.594,
      -0.171, -1.092, -1.854, -2.621, -3.303, -3.801, -3.960, -3.799, -3.151, -2.021, -1.161, -0.834, -0.856, -0.524, -0.885, -1.123, -0.983, 1.450, 1.961, 2.487,
      2.814, 2.452, 2.283, 2.135, 1.659, 1.113, 0.456, -0.148, -1.008, -1.725, -2.395, -3.255, -3.334, -3.509, -3.890, -3.572, -3.054, -2.316, -2.092, -2.151,
      -1.378, 0.000, -0.567, 0.000, 1.619, 1.881, 2.371, 3.079, 2.448, 1.842, 1.195, 0.897, 0.828, 0.356, -0.431, -0.929, -1.451, -2.179, -2.895, -2.751,
      -3.143, -3.354, -3.681, -4.524, -3.337, -2.888, 0.000, 0.000, 0.000, 0.000, 0.000, 0.818, 0.960, 1.182, 1.541, 1.204, 1.035, 0.174, -0.185, -0.493,
      -0.768, -1.307, -1.376, -1.856, -2.268, -2.277, -2.975, -4.125, -3.092, -3.148, -3.235, -3.749, -3.012, -2.000, 0.000, 0.000, 0.000, 0.000, 0.166, -0.379,
      0.018, 0.079, -0.194, -0.486, -1.402, -1.690, -1.363, -2.475, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 1.677, 2.131, 3.213, 3.402, 3.101, 2.915, 2.250, 1.374, 0.290, -0.759, -1.821, -2.654, -3.239, -3.499, -3.835,
      -3.768, -3.587, -2.901, -1.937, -0.988, -1.241, -1.148, -1.206, -1.383, -1.213, -1.581, -1.789, 2.093, 1.882, 2.968, 2.915, 2.810, 2.454, 1.906, 1.122,
      0.148, -0.846, -1.926, -2.748, -3.279, -3.755, -4.090, -3.695, -3.376, -2.701, -2.830, -2.238, -1.685, -2.156, -2.184, -2.097, 0.000, 0.000, 0.000, 1.691,
      1.754, 2.099, 2.907, 2.714, 1.989, 1.398, 0.674, -0.049, -1.293, -1.772, -2.485, -3.091, -3.094, -3.603, -3.437, -3.957, -3.428, -2.453, -2.315, -2.124,
      -1.259, 0.000, 0.000, 0.000, 0.000, 0.000, 1.050, 1.098, 0.959, 1.473, 1.377, 0.913, 0.133, -0.297, -1.162, -1.575, -1.937, -2.740, -2.405, -3.371,
      -3.720, -3.107, -3.614, -3.922, -4.134, -3.449, -2.162, 0.000, -1.181, 0.000, 0.000, 0.000, 0.000, 0.137, -0.243, 0.153, 0.227, -0.243, -0.971, -1.314,
      -1.285, -2.076, -2.846, -1.542, 0.000, -2.383, 0.000, -1.948, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.550, 1.422, 2.899, 3.603, 3.239, 2.594, 1.521, 0.286, -0.995, -2.041, -2.788, -3.212, -3.509, -3.457, -3.071, -3.046, -2.384, -1.601, -1.554, -1.536,
      -1.994, -1.929, -1.814, -1.956, -1.848, -2.053, -2.459, 1.763, 2.485, 3.688, 3.226, 2.888, 2.219, 1.324, 0.124, -1.188, -2.428, -3.151, -3.837, -3.936,
      -3.571, -3.577, -3.310, -2.496, -2.473, -1.423, -1.447, -1.737, -2.686, -2.309, -1.210, 0.000, -0.517, 0.000, 1.927, 1.909, 2.065, 2.800, 2.562, 1.677,
      0.876, -0.066, -1.096, -2.267, -2.780, -3.418, -3.875, -4.097, -3.371, -3.948, -4.360, -2.631, -1.610, -2.868, -1.924, -1.077, 0.000, 0.000, 0.000, 0.000,
      0.000, 1.048, 0.966, 1.095, 1.624, 1.322, 0.579, -0.461, -0.704, -1.626, -2.618, -2.986, -3.171, -4.058, -4.295, -3.122, -3.062, -3.357, -4.425, -3.133,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.095, -0.262, 0.514, 0.462, -0.109, -0.100, -1.341, -1.495, 0.000, -2.125, -3.444, 0.000,
      0.000, 0.000, -2.516, 0.000, 0.000, -1.629, -2.516, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.577, 1.465, 2.711, 3.579, 3.124,
      2.033, 0.649, -0.865, -1.964, -2.671, -3.010, -3.136, -2.683, -2.763, -2.368, -1.971, -1.568, -1.245, -1.238, -1.571, -1.658, -1.942, -1.933, -2.293, -1.652,
      -1.428, -2.932, 1.830, 1.984, 2.859, 3.419, 2.754, 1.834, 0.465, -1.007, -2.190, -3.175, -3.792, -3.697, -3.616, -3.396, -2.469, -1.798, -2.127, -2.512,
      -1.345, -2.457, -1.616, -2.491, -1.170, 0.000, 0.000, 0.000, 0.000, 1.765, 2.133, 2.007, 3.374, 2.400, 1.419, 0.221, -0.893, -2.216, -3.078, -3.872,
      -3.466, -4.024, -3.670, -2.912, -3.174, -2.803, -3.125, 0.000, -2.209, -1.576, 0.000, 0.000, 0.000, 0.000, -0.129, 0.000, 1.230, 0.859, 0.910, 1.701,
      1.311, 0.395, -0.424, -1.738, -2.448, -3.190, -4.320, -3.692, -3.871, -3.410, -3.780, -2.785, -3.792, -3.516, -2.296, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.082, -0.016, 0.534, 0.304, -0.143, -1.236, -0.781, -1.469, -3.212, 0.000, 0.000, -2.263, 0.000, 0.000, 0.000, -2.486, 0.000,
      0.000, 0.000, 0.000, -1.387, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.091, 1.515, 2.945, 3.566, 2.634, 1.156, -0.518, -1.802, -2.509, -2.729,
      -2.657, -2.387, -2.193, -1.790, -1.466, -1.326, -0.859, -0.615, -0.877, -1.012, -1.425, -1.422, -1.286, -1.279, -0.915, -1.493, -1.872, 1.562, 1.884, 3.039,
      3.359, 2.364, 1.092, -0.526, -2.071, -2.915, -3.166, -3.172, -2.748, -2.760, -2.382, -2.163, -1.269, -2.003, -2.002, -2.242, -2.231, -2.412, -1.228, -2.980,
      0.000, 0.000, 0.000, -1.371, 2.071, 2.224, 2.283, 3.039, 2.251, 0.805, -0.741, -2.167, -2.970, -3.493, -3.798, -3.625, -3.794, -3.534, -3.048, -3.354,
      -2.737, -2.761, -1.985, -2.761, 0.000, -2.355, -1.508, -1.102, 0.000, -0.409, 0.000, 1.055, 0.870, 1.204, 1.352, 1.063, 0.211, -1.380, -2.355, -3.161,
      -4.133, -4.764, -3.687, -3.290, -3.761, -3.582, -2.115, -4.045, -3.384, -3.176, -3.009, -2.403, -1.710, -0.611, 0.000, 0.000, 0.000, 0.000, 0.083, 0.099,
      0.451, 0.562, -0.211, -0.838, -1.560, -2.442, -2.401, -4.109, 0.000, 0.000, 0.000, -2.649, -2.649, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.078, 1.600, 3.407, 3.281, 1.694, -0.205, -1.623, -2.291, -2.394, -2.170, -1.784, -1.460, -1.028, -0.817, -0.434,
      -0.301, -0.309, -0.362, -0.469, -0.732, -0.737, -0.579, -1.192, -1.131, -0.837, -1.582, -1.280, 1.980, 1.564, 3.083, 3.118, 1.737, -0.249, -1.808, -2.571,
      -2.726, -2.426, -1.873, -1.678, -1.371, -0.862, -1.120, -1.092, -1.144, -1.802, -1.248, -1.643, -1.835, -1.939, -1.389, -1.859, -2.013, -1.166, 0.000, 1.945,
      2.234, 2.553, 2.903, 1.706, -0.230, -1.950, -2.887, -3.281, -3.347, -2.952, -2.752, -1.994, -2.380, -2.441, -2.706, -2.050, -1.805, -3.030, -1.554, -1.554,
      -1.554, 0.526, 0.000, 0.000, 0.526, 0.000, 1.130, 0.707, 1.142, 1.429, 0.862, -0.787, -2.226, -3.264, -3.620, -3.444, -3.159, -3.688, -3.537, -2.599,
      -2.732, -2.277, -1.410, -3.376, -2.747, -2.460, -2.460, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.127, -0.064, 0.282, 0.017, -0.471, -1.142, -1.947,
      -2.963, -3.179, -2.682, -3.245, -3.625, -2.185, -2.067, -2.290, 0.000, 0.000, -2.067, 0.418, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.604, 1.973, 3.095, 2.182, 0.325, -1.226, -1.851, -1.850, -1.492, -1.138, -0.777, -0.493, -0.189, -0.114, -0.015, -0.088, 0.232, 0.069, 0.043, 0.127,
      -0.304, -0.564, -0.840, -0.470, -0.614, -0.922, -0.658, 2.194, 1.883, 2.872, 2.620, 0.462, -1.437, -2.063, -1.990, -1.604, -1.082, -0.766, -0.375, -0.535,
      -0.512, -0.794, -0.447, -0.580, -0.923, -1.212, -1.735, -1.344, -1.213, -0.657, -1.078, -1.456, -1.233, -2.331, 1.933, 2.025, 2.342, 2.437, 0.474, -1.385,
      -2.343, -2.365, -2.301, -1.718, -1.313, -1.406, -1.053, -0.741, -1.316, -1.029, -1.758, -1.365, -1.365, -1.796, -1.065, -0.267, -0.777, 0.000, 0.000, 0.000,
      0.000, 1.529, 1.168, 1.357, 1.437, 0.004, -1.698, -2.803, -3.021, -2.835, -2.594, -2.237, -2.537, -2.140, -2.258, -2.006, -1.313, -0.215, -1.496, -2.700,
      -2.294, -2.006, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.165, 0.057, 0.174, 0.105, -0.359, -1.756, -2.629, -2.782, -3.136, -2.782, -3.593, -2.494,
      -2.494, 0.000, -1.801, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.104, 1.725, 2.424, 1.037, -0.593,
      -1.379, -1.248, -0.868, -0.584, -0.338, -0.170, 0.040, 0.413, 0.227, 0.446, 0.580, 0.203, -0.013, 0.629, 0.303, 0.723, -0.215, 0.047, 0.527, 1.201,
      0.344, -0.566, 0.971, 1.632, 2.480, 1.270, -0.650, -1.481, -1.183, -0.672, -0.283, 0.219, 0.420, 0.434, 0.479, 0.557, -0.399, -0.398, -0.445, -0.460,
      -0.261, -0.910, -0.781, -1.093, -1.898, -1.966, -0.675, -0.838, -1.592, 3.166, 2.629, 2.139, 1.646, -0.446, -1.651, -1.520, -1.213, -0.894, -0.364, 0.658,
      -0.206, -0.854, -1.195, -0.729, 0.000, 0.000, -0.671, -0.278, -1.164, -1.525, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.569, 1.941, 1.690, 1.038,
      -0.741, -1.870, -1.948, -1.321, -0.834, -0.102, 0.174, -0.268, -1.249, -1.018, -1.212, 0.000, -1.079, 0.000, -2.465, 0.000, 0.000, 0.000, 0.000, -3.851,
      0.000, 0.000, 0.000, 0.410, -0.108, -0.627, 0.112, -1.050, -1.988, -2.384, -2.161, -3.454, -2.943, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.843, 0.255, -0.519, -0.862, -0.212, -1.027, -0.640, 0.000,
      0.000, -0.414, 0.000, -0.991, 0.000, -2.054, 0.000, 0.000, -0.987, 0.000, -1.377, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.312,
      0.107, -0.683, 0.210, 0.000, -0.864, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000,
    };


    std::vector<float> dedx_edges_pl_1 = {
      0.000, 0.500, 1.000, 1.500, 2.000, 2.500, 3.000, 3.500, 4.000, 4.500, 5.000, 5.500, 6.000, 6.500, 7.000, 7.500, 8.000, 9.000, 10.000, 12.000,
      15.000, 20.000, 25.000, 30.000, 35.000, 40.000, 45.000, 50.000,
    };



    std::vector<std::vector<float>> parameters_edges_pl_1 = {
      {0.000, 2.000, 4.000, 7.000, 10.000, 15.000, 20.000, 30.000, 50.000, 100.000, 300.000, 2000.000,},
      {0.000, 0.600, 1.000, 1.500, 3.000, 30.000,}
    };



    std::vector<float> dedx_pdf_pl_1 = {
      1.086, 1.552, 2.578, 2.487, 2.087, 1.732, 1.520, 1.670, 1.483, 1.617, 1.643, 1.662, 1.465, 1.294, 1.237, 1.075, 0.603, -0.012, -0.631, -1.460,
      -2.024, -2.300, -2.142, -2.021, -1.747, -1.316, -0.715, 1.056, 1.971, 2.046, 1.910, 1.950, 1.726, 1.831, 1.531, 1.433, 1.320, 1.288, 0.915, 0.714,
      0.709, 0.473, 0.212, -0.214, -0.651, -1.175, -1.677, -1.945, -1.952, -1.880, -1.723, -1.802, -1.107, -1.092, 1.398, 1.947, 1.624, 1.794, 1.499, 1.609,
      1.279, 0.767, 0.471, 0.858, 0.828, 0.273, -0.039, -0.201, -0.822, -1.096, -0.910, -1.281, -1.452, -2.103, -2.170, -2.191, -1.847, -1.980, -2.008, -2.268,
      -1.980, 1.024, 0.829, 0.605, 0.846, 0.794, 0.416, 0.206, -0.202, 0.105, -0.696, -0.409, -0.802, -1.262, -1.957, -1.957, -2.191, -2.176, -1.994, -2.347,
      -3.165, -3.153, -3.616, -3.092, 0.000, 0.000, -1.013, 0.000, 0.159, -0.152, 0.066, -0.071, -0.290, -1.126, -1.693, -0.631, -1.627, -2.449, -1.661, -1.843,
      -1.438, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.325, 2.765, 3.073, 3.016, 2.311,
      2.175, 2.110, 2.167, 1.905, 1.778, 1.384, 1.007, 0.412, -0.144, -0.835, -1.315, -2.158, -2.918, -3.249, -3.440, -2.881, -1.816, -1.361, -1.174, -0.921,
      -1.204, -1.789, 1.410, 2.913, 2.785, 2.541, 2.240, 2.093, 1.942, 1.911, 1.736, 1.273, 0.916, 0.585, 0.293, -0.569, -0.808, -1.220, -1.635, -2.336,
      -2.877, -2.868, -2.534, -2.046, -1.632, -1.836, -1.927, -1.836, -0.989, 1.198, 2.234, 1.642, 2.311, 2.302, 2.161, 1.390, 1.187, 0.853, 0.593, 0.328,
      -0.015, -0.380, -0.737, -0.940, -1.397, -2.101, -2.158, -2.265, -2.915, -2.529, -2.017, -2.082, -2.964, -1.389, -0.896, -0.513, 0.860, 1.020, 0.898, 0.984,
      1.202, 0.873, 0.679, 0.377, -0.019, -0.082, -1.063, -1.014, -1.270, -2.113, -2.567, -1.871, -2.053, -3.375, -3.231, -3.291, -3.068, -2.931, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.178, -0.077, -0.126, 0.221, -0.686, -0.771, -0.925, -2.450, -1.442, -1.603, -1.351, -1.351, -1.970, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.552, 2.618, 3.403, 2.939, 2.673, 2.550, 2.325, 1.926, 1.357, 0.712,
      -0.014, -0.806, -1.586, -2.187, -2.786, -3.099, -3.479, -3.601, -3.408, -2.824, -2.276, -2.102, -1.944, -1.463, -1.189, -1.338, -0.827, 1.448, 2.641, 2.919,
      3.178, 2.547, 2.446, 2.023, 1.803, 1.315, 0.494, -0.155, -0.883, -1.556, -2.022, -2.332, -3.026, -3.355, -3.256, -3.293, -3.444, -2.261, -2.277, -1.793,
      -1.820, -1.515, 0.000, 0.026, 1.394, 2.147, 2.204, 2.230, 2.192, 1.943, 1.477, 1.256, 0.763, 0.276, -0.239, -1.044, -1.289, -1.799, -2.130, -3.027,
      -3.107, -3.222, -2.875, -3.081, -2.304, -1.610, -0.929, -1.016, 0.000, 0.000, 0.000, 0.778, 0.856, 0.971, 1.269, 1.161, 1.138, 0.500, -0.253, -0.322,
      -0.931, -1.491, -1.987, -2.320, -2.669, -2.680, -2.511, -2.360, -4.364, -4.221, -2.972, 0.000, -3.059, 0.000, 0.000, 0.000, 0.000, 0.000, 0.113, -0.208,
      0.306, 0.480, -0.057, -0.032, -0.703, -1.288, -1.766, -1.697, 0.000, -1.766, 0.000, 0.000, 0.000, 0.000, -2.599, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 1.262, 2.356, 3.375, 3.187, 2.870, 2.576, 1.994, 1.235, 0.315, -0.635, -1.409, -2.302, -2.932, -3.172, -3.468,
      -3.408, -3.351, -3.013, -2.873, -2.169, -1.747, -2.219, -1.698, -1.174, -0.290, -0.610, -2.052, 1.984, 2.107, 2.838, 3.356, 2.614, 2.463, 1.996, 1.157,
      0.421, -0.683, -1.618, -2.346, -3.021, -3.268, -3.816, -4.299, -3.863, -3.371, -3.285, -2.537, -1.984, -1.026, -1.923, -1.464, 0.000, -0.770, 0.000, 1.598,
      2.047, 2.449, 2.533, 2.348, 2.144, 1.520, 0.847, 0.002, -0.910, -1.349, -2.224, -2.561, -3.908, -3.184, -3.362, -3.465, -3.760, -2.771, -1.921, -1.126,
      -0.944, -0.789, 0.000, 0.000, 0.000, 0.000, 0.841, 0.934, 0.742, 1.463, 1.160, 0.780, 0.394, -0.512, -0.807, -1.360, -1.948, -3.537, -2.785, -3.531,
      -3.425, -4.604, -3.320, -3.594, -3.468, -4.217, 0.000, 0.000, -1.522, 0.000, 0.000, 0.000, 0.000, 0.101, -0.238, 0.413, 0.401, 0.092, -0.585, -1.096,
      -0.907, -2.656, -1.772, -1.823, -0.842, 0.000, 0.000, 0.000, 0.000, -2.612, 0.000, 0.000, 0.000, 0.000, -0.214, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.941, 2.256, 3.061, 3.237, 2.853, 2.245, 1.314, 0.271, -0.870, -1.671, -2.427, -2.905, -3.055, -3.266, -3.146, -2.916, -2.826, -2.515, -1.942, -1.710,
      -1.628, -1.557, -1.221, -1.196, -0.953, -1.036, -1.392, 1.886, 2.216, 2.938, 3.360, 2.803, 2.268, 1.322, 0.277, -0.765, -1.859, -2.710, -3.315, -3.418,
      -3.510, -3.266, -3.446, -3.580, -2.637, -2.451, -2.018, -1.952, -1.721, -1.421, -1.604, -0.911, 0.000, 1.287, 1.421, 2.089, 1.851, 2.934, 2.600, 1.713,
      1.043, 0.052, -0.850, -1.996, -2.540, -3.235, -3.209, -3.604, -3.501, -4.175, -4.237, -3.177, -3.354, -1.656, -2.421, -1.492, -0.511, 0.000, 0.000, 0.000,
      0.000, 0.797, 0.786, 0.817, 1.432, 1.126, 0.623, 0.069, -0.627, -1.450, -2.002, -2.460, -3.766, -3.410, -3.465, -3.533, -3.856, -2.898, -3.578, -3.039,
      -2.768, -1.730, 0.000, 0.000, 0.005, 0.000, 0.005, 0.000, 0.086, -0.159, 0.611, 0.161, -0.042, -0.173, -0.451, -0.989, -3.186, -1.745, 0.000, -2.834,
      -1.773, -2.794, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.030, 1.923, 3.060, 3.195, 2.695,
      1.745, 0.519, -0.691, -1.620, -2.314, -2.630, -2.750, -2.746, -2.723, -2.558, -2.527, -2.100, -1.908, -1.690, -1.481, -1.504, -1.832, -1.652, -3.179, -2.403,
      -2.358, -1.752, 2.332, 1.853, 2.813, 3.261, 2.609, 1.822, 0.579, -0.658, -1.845, -2.669, -3.178, -3.427, -3.535, -2.830, -2.818, -2.466, -2.976, -3.036,
      -2.837, -1.892, -1.862, -1.847, 0.000, -1.362, -0.668, 0.000, 0.025, 1.734, 2.189, 2.361, 3.013, 2.439, 1.455, 0.337, -1.041, -1.927, -3.095, -3.508,
      -3.799, -3.685, -3.210, -2.957, -2.397, -2.333, -2.502, -2.613, -3.293, -1.374, -1.173, -1.684, 0.000, 0.000, 0.000, 0.000, 0.811, 0.579, 0.856, 1.443,
      1.254, 0.423, -0.225, -1.360, -2.023, -2.783, -3.525, -2.908, -4.457, -3.659, -3.625, -4.348, -4.225, -2.416, -2.556, -2.753, -2.416, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.109, -0.191, 0.517, 0.660, 0.130, -0.296, -0.771, -1.767, -2.104, -2.461, -2.461, 0.000, 0.000, 0.000, 0.000, 0.000, -2.317,
      0.000, 0.000, 0.000, 0.248, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.409, 2.414, 3.287, 3.049, 2.152, 0.899, -0.426, -1.440, -2.152, -2.465,
      -2.542, -2.470, -2.428, -2.104, -1.854, -1.876, -1.554, -1.265, -1.320, -1.236, -0.991, -1.209, -2.589, -0.816, -0.816, -1.902, -1.240, 1.732, 2.014, 3.164,
      3.128, 2.337, 1.113, -0.371, -1.632, -2.561, -2.878, -3.077, -2.662, -2.467, -2.390, -2.412, -2.442, -1.848, -2.209, -2.278, -2.048, -1.842, -1.231, -1.263,
      -0.926, -0.638, -0.926, -0.233, 1.858, 2.037, 2.104, 2.607, 2.248, 0.985, -0.533, -1.749, -2.736, -3.122, -3.459, -3.106, -2.894, -2.893, -2.768, -2.837,
      -2.613, -2.646, -2.950, -1.815, -1.564, 0.000, -1.920, -1.409, 0.000, 0.000, 0.000, 0.772, 0.622, 1.025, 1.518, 1.065, 0.190, -0.963, -2.071, -2.746,
      -3.030, -2.941, -3.180, -3.412, -4.256, -3.058, -3.136, -2.005, -2.664, -3.229, 0.000, -1.648, -0.396, 0.000, 0.000, 0.000, 0.000, 0.000, 0.140, -0.077,
      0.259, 0.211, -0.097, -0.962, -1.470, -2.539, -2.869, -1.968, -3.734, -3.623, -2.622, -2.930, -2.495, -2.043, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 1.067, 2.097, 3.066, 2.582, 1.301, -0.113, -1.229, -1.893, -2.116, -2.106, -1.943, -1.726, -1.487, -1.291, -1.231,
      -1.128, -1.025, -0.778, -0.904, -0.663, -0.860, -1.005, -1.151, -1.295, -1.475, -1.247, -0.261, 2.101, 2.301, 3.100, 2.777, 1.586, -0.031, -1.419, -2.234,
      -2.550, -2.447, -2.076, -1.820, -1.634, -1.458, -1.386, -1.217, -1.436, -1.365, -1.806, -1.383, -1.566, -1.274, -1.839, -0.997, -2.270, -1.354, -2.452, 1.739,
      1.741, 2.195, 2.655, 1.464, -0.112, -1.579, -2.484, -2.761, -2.817, -2.770, -2.201, -2.629, -2.200, -1.897, -2.808, -1.505, -1.012, -2.121, -2.218, -2.275,
      -1.058, 0.000, 0.000, 0.000, 0.000, -0.771, 1.130, 0.961, 1.101, 1.403, 0.760, -0.621, -1.909, -2.638, -3.238, -3.495, -3.468, -3.664, -3.099, -2.762,
      -2.150, -3.202, -1.949, -1.793, -1.793, -1.698, -1.921, -1.698, 0.000, 0.000, 0.000, 0.000, 0.000, 0.139, -0.023, 0.317, 0.279, -0.161, -1.450, -1.972,
      -2.342, -3.319, -2.938, -2.285, 0.000, -2.808, -1.170, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.611, 0.000, 0.000, 0.000, 0.000,
      0.683, 2.224, 2.594, 1.644, 0.265, -0.872, -1.502, -1.669, -1.572, -1.293, -1.217, -0.881, -0.668, -0.454, -0.532, -0.271, -0.213, -0.108, 0.069, 0.135,
      -0.094, -0.270, -0.452, -0.052, 0.032, -0.607, -0.561, 2.349, 2.162, 2.803, 2.153, 0.447, -1.068, -1.845, -1.897, -1.717, -1.402, -1.164, -0.903, -0.796,
      -0.775, -0.540, -1.094, -0.610, -0.557, -0.564, -0.998, -0.915, -0.866, -1.020, 0.471, -0.104, -0.355, -0.222, 1.888, 2.024, 2.175, 2.190, 0.352, -1.218,
      -2.022, -2.134, -2.249, -1.892, -1.674, -1.533, -1.778, -1.650, -1.104, -1.450, -0.870, -1.209, -1.555, -1.245, -3.355, -2.749, -0.957, 0.000, 0.000, -1.650,
      0.000, 1.356, 1.218, 1.125, 1.138, 0.156, -1.475, -2.387, -2.718, -2.966, -2.990, -2.042, -2.340, -2.944, -1.868, -2.656, -1.935, -3.080, -3.167, -2.944,
      0.011, -1.781, -2.068, -0.970, 0.000, -2.068, 0.000, 0.000, 0.154, -0.008, 0.294, 0.206, -0.646, -1.361, -2.015, -2.746, -3.170, -2.993, -3.198, -3.304,
      -3.080, 0.000, -2.387, 0.000, -2.387, 0.000, -1.289, 0.000, -0.595, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.192, 1.814, 1.554, 0.716, -0.368,
      -1.040, -1.173, -0.999, -0.748, -0.485, -0.387, -0.212, -0.265, 0.806, 0.164, 0.639, 0.410, 0.780, 0.455, 0.658, 0.335, 0.084, -0.337, -0.421, -0.818,
      -0.512, -0.895, 3.264, 1.763, 2.385, 0.971, -0.474, -1.226, -1.160, -0.830, -0.560, -0.428, -0.427, -0.027, 0.232, 0.087, 0.240, 0.896, 0.764, 1.233,
      1.249, 0.172, 1.354, 0.415, 0.910, 0.000, 0.000, -1.281, -0.856, 2.007, 2.817, 2.316, 1.047, -0.599, -1.254, -1.247, -1.056, -0.524, -0.847, -0.107,
      -0.618, -0.841, -0.457, -0.359, 0.136, -1.144, 0.060, -0.892, -1.242, -1.180, -2.191, 0.000, 0.000, 0.000, 0.000, 0.000, 1.796, 1.295, 0.555, 0.654,
      -0.558, -1.364, -1.697, -1.757, -1.167, -1.232, -1.330, -0.933, -1.145, -2.144, -0.912, -1.199, 0.000, -2.144, -0.871, 0.000, -1.787, -3.397, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.180, 0.351, -0.368, 0.418, -0.548, -1.438, -1.479, -1.485, -2.728, -2.140, 0.000, -2.392, -2.546, 0.000, 0.000, -3.239, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.750, 0.939, 0.030, -0.500, -1.426, -1.403, -1.612, -1.805,
      -1.645, -1.415, -0.215, 0.000, -1.361, -0.841, -1.645, 0.000, -0.527, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.774, 1.019,
      0.187, -0.513, -0.630, -0.190, -0.690, 0.000, 0.000, -0.260, 0.000, -0.881, 0.000, 0.000, -1.443, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -2.510, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -7.746, 0.000, 0.000, 0.000, 0.000, -2.191, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.343, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000,
    };


    std::vector<float> dedx_edges_pl_2 = {
      0.000, 0.500, 1.000, 1.500, 2.000, 2.500, 3.000, 3.500, 4.000, 4.500, 5.000, 5.500, 6.000, 6.500, 7.000, 7.500, 8.000, 9.000, 10.000, 12.000,
      15.000, 20.000, 25.000, 30.000, 35.000, 40.000, 45.000, 50.000,
    };



    std::vector<std::vector<float>> parameters_edges_pl_2 = {
      {0.000, 2.000, 4.000, 7.000, 10.000, 15.000, 20.000, 30.000, 50.000, 100.000, 300.000, 2000.000,},
      {0.000, 0.600, 1.000, 1.500, 3.000, 30.000,}
    };



    std::vector<float> dedx_pdf_pl_2 = {
      0.379, 1.402, 2.485, 2.947, 2.220, 1.720, 1.311, 1.055, 1.142, 1.403, 1.672, 1.899, 1.963, 1.708, 1.604, 1.571, 1.285, 1.002, 0.461, -0.932,
      -2.118, -2.430, -2.494, -2.335, -1.742, -1.146, -0.893, 0.717, 2.388, 2.302, 2.479, 3.131, 2.113, 2.728, 2.442, 2.332, 2.297, 1.688, 1.871, 1.360,
      1.412, 1.050, 1.058, 0.600, 0.140, -0.425, -1.266, -2.031, -2.053, -1.884, -1.741, -1.343, -1.083, -0.555, 0.816, 2.286, 1.785, 1.838, 1.967, 1.816,
      2.350, 1.991, 1.587, 1.716, 0.677, 1.030, 0.320, 0.463, 0.463, 0.303, 0.034, -0.654, -1.049, -1.883, -2.466, -1.898, -2.091, -3.631, -1.979, -0.789,
      -2.340, 0.781, 1.125, 0.915, 0.744, 1.125, 0.932, 0.587, 0.434, 0.588, 0.091, -0.341, -0.550, -0.459, -1.139, -0.867, -3.237, -1.427, -2.526, -2.508,
      -2.402, -3.943, -3.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.073, -0.014, 0.021, 0.107, 0.165, -1.339, -1.222, -1.157, 0.000, -1.445, -1.781, -0.529,
      0.000, 0.000, 0.165, 0.000, 0.000, -0.934, 0.000, 0.165, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.677, 3.761, 3.114, 2.252,
      1.702, 1.608, 1.988, 2.438, 2.326, 2.054, 1.624, 1.042, 0.494, 0.060, -0.441, -1.274, -2.575, -3.822, -4.104, -2.932, -1.467, -0.778, -0.383, 0.055,
      -0.443, -0.232, 1.134, 2.190, 2.875, 2.733, 3.210, 3.295, 2.555, 3.219, 3.104, 2.665, 2.446, 1.828, 1.339, 1.190, 0.589, 0.333, -0.799, -1.612,
      -2.782, -3.133, -2.841, -3.304, -2.023, -1.168, -0.540, -1.350, -1.350, 1.281, 1.953, 2.169, 0.000, 3.469, 3.001, 2.960, 3.221, 2.782, 1.949, 1.151,
      1.407, 0.578, 0.033, 0.183, -0.685, -1.141, -1.449, -3.053, -3.344, -2.593, -3.024, -2.513, -2.330, -1.974, 0.000, 0.000, 0.800, 0.895, 0.887, 1.502,
      1.789, 1.461, 1.183, 0.850, 0.665, 0.691, -0.002, -0.139, 0.180, -1.073, -1.612, -1.717, -2.218, -2.045, -3.422, -3.663, -3.059, -2.911, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.027, -0.132, 0.629, 0.572, 0.738, -0.506, 0.070, -0.980, -0.064, -0.757, -2.010, 0.000, -0.757, -1.673, 0.000, -0.064, 0.000,
      -1.450, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.403, 4.215, 2.994, 2.341, 2.318, 2.670, 2.590, 1.979, 1.260,
      0.538, -0.247, -0.848, -1.662, -2.580, -3.616, -4.429, -4.736, -4.283, -2.430, -1.280, -0.585, -0.433, -0.541, -0.305, 0.101, 0.389, 0.449, 2.569, 2.396,
      3.383, 3.703, 4.168, 4.064, 3.494, 2.691, 2.231, 1.310, 0.563, -0.241, -0.957, -1.923, -2.457, -3.441, -4.523, -4.093, -4.288, -3.843, -2.067, 0.000,
      0.000, -0.362, 0.331, 0.000, 1.779, 1.371, 2.798, 4.407, 3.586, 3.753, 3.758, 3.076, 2.683, 1.411, 1.143, 0.174, -0.380, -1.459, -2.265, -2.231,
      -2.837, -3.324, -3.899, -3.417, -2.655, 0.000, 0.000, -0.603, 0.000, 0.000, 0.000, 0.914, 1.110, 0.725, 1.410, 2.387, 1.170, 1.438, 0.910, 0.592,
      -0.160, -0.647, -1.860, -0.986, -2.665, -3.217, -2.300, -3.235, -4.058, -4.443, -4.913, -3.997, 0.000, -2.428, 0.000, 0.000, 0.000, 0.000, 0.018, -0.118,
      0.875, 0.251, 0.394, 0.251, -0.008, -0.259, 0.434, -1.176, -1.646, 0.000, -1.646, -1.646, -0.547, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, -0.840, 3.254, 3.505, 3.091, 2.740, 3.073, 2.727, 1.690, 0.815, 0.009, -0.873, -1.741, -2.805, -3.763, -4.446,
      -4.612, -4.406, -3.288, -1.937, -1.451, -1.225, -1.063, -1.218, -0.495, -0.453, 0.189, -0.042, 1.376, 1.686, 1.676, 2.705, 4.347, 4.510, 4.027, 2.928,
      2.067, 0.821, -0.365, -1.681, -2.780, -3.980, -4.523, -6.732, -4.936, -4.453, -3.491, -3.871, 0.000, -2.130, 0.000, 0.000, -0.521, 0.000, 0.000, 1.372,
      1.596, 2.762, 3.801, 4.389, 3.747, 2.872, 1.788, 1.347, 1.028, -0.373, -1.452, -2.200, -2.917, -4.142, -4.187, -4.292, -3.775, -3.264, -2.771, 0.000,
      0.000, -1.067, 0.000, 0.000, 0.000, 0.000, 0.938, 0.703, 0.800, 1.619, 1.845, 1.425, 0.994, -0.055, -0.585, 0.072, -1.943, -2.088, -3.001, -3.768,
      -3.008, -3.768, -4.140, -3.592, -3.362, -3.425, 0.000, -1.027, 0.000, 0.000, 0.000, 0.000, 0.000, 0.011, 0.026, 0.766, 0.658, 0.494, -0.237, 0.138,
      -0.555, -1.066, -1.654, -1.249, -2.347, 0.000, 0.000, -0.555, 0.000, -2.347, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      -0.246, 2.152, 3.095, 3.257, 3.264, 3.106, 1.859, 0.724, -0.209, -1.222, -2.461, -3.482, -3.992, -4.069, -3.763, -3.125, -2.295, -1.587, -1.529, -1.096,
      -1.535, -1.946, -2.180, -2.811, -1.321, -1.162, -0.988, 1.739, 1.419, 2.332, 2.933, 4.018, 4.044, 3.495, 1.873, 0.574, -0.490, -2.327, -3.995, -4.308,
      -5.123, -5.138, -5.047, -3.675, 0.000, -2.749, -2.509, -3.518, -2.537, 0.000, 0.000, 0.000, 0.000, 0.000, 1.139, 1.022, 1.745, 2.048, 3.068, 3.723,
      2.669, 1.376, 0.192, -1.032, -1.879, -3.306, -3.494, -3.571, -3.690, -3.460, -4.889, -3.335, -3.192, -3.622, -3.111, 0.000, 0.000, -2.418, 0.000, 0.000,
      0.000, 1.283, 0.993, 1.050, 2.504, 1.712, 1.567, 0.555, -0.523, -0.571, -1.560, -2.795, -3.931, -3.854, -4.037, -3.868, 0.000, -3.944, -3.275, 0.000,
      0.000, -2.245, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.014, 0.086, 0.682, 0.770, -0.190, -0.615, 1.063, -1.289, -1.576, -0.190, -1.576, 0.000,
      0.000, -2.493, 0.000, 0.000, -1.576, 0.000, -1.982, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.182, 0.919, 3.219, 3.521, 3.576,
      2.518, 0.907, -0.221, -1.448, -2.803, -3.636, -3.800, -3.423, -2.841, -2.399, -2.021, -1.425, -0.904, -0.944, -1.420, -1.884, -1.982, -1.674, -2.836, -1.977,
      -2.878, -2.696, 1.338, 1.905, 1.844, 3.334, 4.429, 4.117, 2.278, 0.668, -0.769, -2.746, -3.895, -4.581, -4.745, -4.510, -3.451, -3.991, -3.803, 0.000,
      -3.216, -3.611, 0.000, -1.270, 0.000, 0.000, 0.000, 0.000, 0.000, 1.263, 1.831, 1.656, 1.892, 3.636, 3.556, 2.347, 0.563, -1.207, -2.300, -3.702,
      -3.550, -4.059, 0.000, -4.438, -2.551, -3.832, 0.000, -2.733, -2.887, 0.000, 0.000, -0.941, 0.000, 0.000, 0.000, 0.000, 1.128, 0.585, 1.090, 2.955,
      1.530, 1.668, 0.380, -0.815, -1.866, -3.043, -4.356, -3.280, -2.952, -3.257, 0.000, -2.998, 0.000, -2.798, -3.309, 0.000, -2.392, -1.699, 0.000, 0.000,
      0.000, 0.000, 0.000, -0.004, 0.208, 1.091, 0.537, -0.027, 0.189, 0.412, -0.974, -2.113, -2.807, 0.000, -2.989, 0.000, 0.000, -1.890, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.155, 1.461, 3.355, 4.006, 3.426, 1.581, -0.079, -1.536, -2.843, -3.450,
      -3.248, -2.858, -2.349, -1.827, -1.517, -0.948, -0.562, -0.444, -0.494, -0.868, -1.011, -1.038, -1.360, -1.366, -1.154, -0.895, -2.005, 2.193, 1.346, 1.751,
      3.397, 3.958, 3.115, 1.424, -0.577, -2.590, -3.547, -3.795, -3.569, -3.427, -3.211, -2.610, -2.213, -2.486, -2.538, -2.962, -2.538, -1.893, -2.104, 0.000,
      -2.643, 0.000, 0.000, 0.000, 1.544, 0.904, 0.884, 3.102, 4.381, 3.402, 1.646, -0.562, -2.582, -3.862, -4.419, -4.434, -3.792, -2.875, -2.875, 0.000,
      -2.182, -1.084, -2.588, 0.000, 0.000, -1.489, 0.000, 0.000, 0.000, 0.000, 0.000, 0.998, 0.527, 0.834, 2.198, 2.209, 1.245, 0.589, -1.658, -2.895,
      -3.525, -4.513, -4.172, -4.107, -4.618, -3.009, -3.925, -4.262, 0.000, -2.721, -3.009, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.005, -0.023,
      2.052, 0.484, 0.314, -0.263, -0.426, -3.354, -3.521, -3.441, -1.937, 0.000, -0.956, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, -0.633, 1.735, 3.516, 4.141, 2.574, 0.158, -1.687, -2.822, -2.918, -2.510, -2.084, -1.688, -1.336, -0.871, -0.552,
      -0.241, -0.147, -0.152, -0.228, -0.365, -0.402, -0.546, -1.028, -0.854, -1.554, -0.987, -1.715, 0.517, 1.561, 1.175, 4.022, 4.018, 2.024, -0.491, -2.552,
      -3.238, -3.017, -3.002, -2.311, -1.596, -1.939, -1.186, -1.556, -1.757, -0.924, -1.287, -1.657, -2.555, -1.862, -2.738, 0.000, -1.757, -1.757, -1.064, 1.118,
      1.183, 1.457, 2.701, 2.885, 1.449, -0.703, -2.532, -3.321, -3.464, -2.741, -2.698, -3.219, -2.015, -0.917, -2.932, -0.734, -2.526, -1.833, -2.526, -2.526,
      -1.428, 0.000, 0.000, 0.000, 0.000, 0.000, 0.900, 0.759, 1.450, 2.159, 1.434, 0.585, -1.239, -2.557, -3.728, -3.983, -4.249, -3.802, -3.109, -4.585,
      -3.332, -3.738, -4.719, -1.946, -4.431, 0.000, -2.234, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.035, 1.189, 0.270, 2.413, -0.196, -1.041, -1.011,
      -1.205, -1.829, -1.647, 0.000, 0.000, 0.000, -2.746, 0.000, 0.000, 0.000, -2.746, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      -0.071, 2.525, 3.407, 3.450, 0.781, -1.650, -2.463, -2.194, -1.703, -1.378, -1.021, -0.715, -0.453, -0.293, -0.034, 0.017, 0.219, 0.102, 0.216, 0.206,
      0.012, -0.188, -0.580, -0.553, -0.967, -0.869, -1.057, 1.658, 0.160, 0.144, 3.752, 2.645, -0.160, -2.332, -2.543, -2.318, -1.414, -0.956, -1.448, -0.524,
      -0.260, 0.677, -0.052, 0.000, 0.000, -0.658, -0.636, -1.191, -2.267, -1.980, -0.801, 0.000, -1.420, -2.267, 1.383, -0.211, 1.966, 2.760, 0.967, 0.011,
      -2.010, -2.844, -2.859, -1.977, -1.821, -0.912, -1.856, -1.722, -2.703, 0.000, -2.010, 0.000, 0.000, -1.787, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.431, -0.350, 0.541, 0.835, 1.992, -1.179, -2.137, -2.012, -2.487, -2.430, 0.000, 0.000, 0.000, 0.000, 0.000, -5.321, 0.000, 0.000, -4.222,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.009, -0.343, 0.000, 1.345, 0.580, -0.265, -2.022, 0.000, -2.022, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.470, 1.583, 2.424, 1.682, -0.980,
      -1.937, -1.384, -1.002, -0.708, -0.384, -0.182, 0.016, -0.032, 0.154, 0.202, 0.492, 0.350, 0.268, 0.443, 0.574, 0.561, 0.671, 0.138, -0.342, -0.411,
      -0.251, -0.411, 0.000, -2.011, -2.694, 2.108, 2.953, -0.236, -0.711, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.513, -2.365, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.485, -1.149, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, -0.091, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.660, 0.357, -0.949, -1.272, -0.952, -0.794, -0.486, -0.258,
      0.287, 0.481, -0.463, 0.678, -0.187, 0.000, 0.000, 0.376, 0.715, 0.000, 0.000, -0.145, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
      0.000, 0.000, 0.000, 0.000, 0.000,
    };
  };
}


#endif
