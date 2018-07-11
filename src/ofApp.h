#pragma once

#include "ofMain.h"
#include "FluidSolver.hpp"
#include "ofxGui.h"

class ofApp : public ofBaseApp
{
public:
    void setup();
    void update();
    void draw();
    void keyPressed(int key);
    void mouseDragged(int x, int y, int button);

private:
    FluidSolver solver;
    
    ofxPanel gui;
    ofParameter<bool> drawVelocity;
    ofParameter<bool> drawDensity;
    ofParameter<bool> addForce;
    ofParameter<bool> addSource;
};
