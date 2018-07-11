#pragma once

#include <iostream>
#include <vector>
#include <array>

#include "ofMain.h"

class FluidSolver
{
public:
    void setup(int w, int h, int winW, int winH);
    void update();
    void drawVelocity();
    void drawDensity();
    void addForce(int x, int y);
    void addSource(int x, int y);
    int getIndexOfCell(int x, int y);

private:
    void addSource(std::vector<float>& x, const std::vector<float>& s);
    void setBound(int bound, std::vector<float>& x);
    void linearSolve(int bound, std::vector<float>& x, std::vector<float>& x0, float a, float c);
    void diffuse(int bound, std::vector<float>& x, std::vector<float>& x0, float diff);
    void advect(int bound, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& u, std::vector<float>& v);
    void project(std::vector<float>& u, std::vector<float>& v, std::vector<float>& p, std::vector<float>& div);
    void velocityStep();
    void densityStep();

    // Utils
    inline void swap(std::vector<float>& x0, std::vector<float>& x);
    inline int idx(int i, int j);
    inline ofVec2f win2Cell(int x, int y);

    int winW, winH;
    int cellW, cellH;
    int numCells;

    float dt;
    float diff, viscosity;
    float force, source;
    
    int numItr;
    
    std::vector<float> density, densityPrev;
    std::vector<float> u, v, uPrev, vPrev;
    
    ofVec2f mousePrev;
};
