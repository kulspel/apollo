function [NodesAB,NodesBC,NodesCD,NodesDA,NodesTubes] = BCNodes()
NodesAB = [3,   9, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, ...
          255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, ...
          271, 272, 273, 274, 275, 276, 277, 278, 279]

NodesBC = [7,   9, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, ...
          136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, ...
          152, 153];

NodesCD = [2,   7, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117];

NodesDA = [2,  3, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, ...
          31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, ...
          47, 48];

NodesTubes = [10,  11,  12, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, ...
              199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 280, 281, 282, 283, 284, ...
              285, 286, 287, 288, 289, 290, 291];
end