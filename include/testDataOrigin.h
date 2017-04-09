//
//
//

#ifndef INITIAL_DATA1_H
#define INITIAL_DATA1_H

// Transform refPoints to create point set:
// translateX: 0.0
// translateY: 0.0
// translateZ: 0.0
//
// rotateX: 0.0
// rotateY: 0.0
// rotateZ: 0.0
// rotateOrder: xyz
//
// scaleX: 1.0
// scaleY: 1.0
// scaleZ: 1.0


// Number of samples
static
const int testDataOrigin1_num = 26;


// Values
static
double testDataOrigin1_values[26*3] = {
    -0.5, -0.5, 0.5,
    0.0, -0.5, 0.5,
    0.5, -0.5, 0.5,
    -0.5, 0.0, 0.5,
    0.0, 0.0, 0.5,
    0.5, 0.0, 0.5,
    -0.5, 0.5, 0.5,
    0.0, 0.5, 0.5,
    0.5, 0.5, 0.5,
    -0.5, 0.5, 0.0,
    0.0, 0.5, 0.0,
    0.5, 0.5, 0.0,
    -0.5, 0.5, -0.5,
    0.0, 0.5, -0.5,
    0.5, 0.5, -0.5,
    -0.5, 0.0, -0.5,
    0.0, 0.0, -0.5,
    0.5, 0.0, -0.5,
    -0.5, -0.5, -0.5,
    0.0, -0.5, -0.5,
    0.5, -0.5, -0.5,
    -0.5, -0.5, 0.0,
    0.0, -0.5, 0.0,
    0.5, -0.5, 0.0,
    0.5, 0.0, 0.0,
    -0.5, 0.0, 0.0
};


#endif // INITIAL_DATA1_H
