#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup()
{
    ofSetBackgroundColor(0);
    
    solver.setup(100, 100, ofGetWidth(), ofGetHeight());
    
    gui.setup();
    gui.add(drawVelocity.set("drawVelocity", true));
    gui.add(drawDensity.set("drawDensity", false));
    gui.add(addForce.set("addForce", true));
    gui.add(addSource.set("addSource", false));
}

//--------------------------------------------------------------
void ofApp::update()
{
    ofSetWindowTitle(ofToString(ofGetFrameRate()));
    
    solver.update();
}

//--------------------------------------------------------------
void ofApp::draw()
{
    if(drawVelocity)
        solver.drawVelocity();
    else if(drawDensity)
        solver.drawDensity();
    
    gui.draw();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key)
{
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button)
{
    if(addForce)
        solver.addForce(x, y);
    else if(addSource)
        solver.addSource(x, y);
}

//--------------------------------------------------------------
