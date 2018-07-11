#include "FluidSolver.hpp"

//--------------------------------------------------------------
// Public
//--------------------------------------------------------------

void FluidSolver::setup(int w, int h, int winW, int winH)
{
    this->cellW = w - 2;
    this->cellH = h - 2;
    
    this->winW = winW;
    this->winH = winH;
    
    numCells = (cellW + 2) * (cellH + 2);
    
    dt = 0.1;
    diff = 0.0;
    viscosity = 0.0;
    
    force = 5.0;
    source = 100.0;
    
    numItr = 20;
    
    density = std::vector<float>(numCells);
    densityPrev = std::vector<float>(numCells);
    u = std::vector<float>(numCells);
    uPrev = std::vector<float>(numCells);
    v = std::vector<float>(numCells);
    vPrev = std::vector<float>(numCells);
}

void FluidSolver::update()
{
    for(int i = 0; i < numCells; i++)
    {
        uPrev[i] = vPrev[i] = densityPrev[i] = 0.0;
    }
    
    velocityStep();
    densityStep();
}

void FluidSolver::drawVelocity()
{
    ofPushMatrix();
    ofVec2f aspect = ofVec2f(float(winW) / (float(cellW) + 2), float(winH) / (float(cellH) + 2));
    ofScale(aspect.x, aspect.y);
    for(int i = 0; i <= cellW; i++)
    {
        float x = (i + 1.0);
        for(int j = 0; j <= cellH; j++)
        {
            float y = (j + 1.0);
            ofDrawLine(x, y, x + u[idx(i, j)], y + v[idx(i, j)]);
        }
    }
    ofPopMatrix();
}

void FluidSolver::drawDensity()
{
    ofPushMatrix();
    ofVec2f aspect = ofVec2f(float(winW) / float((cellW + 2)), float(winH) / (float(cellH + 2)));
    ofScale(aspect.x, aspect.y);

    glBegin(GL_QUADS);
    for(int i = 0; i <= cellW; i++)
    {
        float x = i + 0.5;
        for(int j = 0; j <= cellH; j++)
        {
            float y = j + 0.5;

            float d00 = density[idx(i,   j)];
            float d01 = density[idx(i,   j+1)];
            float d10 = density[idx(i+1, j)];
            float d11 = density[idx(i+1, j+1)];
            
            glColor3f(d00, d00, d00); glVertex2f(x, y);
            glColor3f(d10, d10, d10); glVertex2f(x+1, y);
            glColor3f(d11, d11, d11); glVertex2f(x+1, y+1);
            glColor3f(d01, d01, d01); glVertex2f(x, y+1);
        }
    }
    glEnd();
    ofPopMatrix();
}

void FluidSolver::addForce(int x, int y)
{
    int index = getIndexOfCell(x, y);
    u[index] = force * (x - mousePrev.x);
    v[index] = force * (mousePrev.y - y);
    
    mousePrev = ofVec2f(x, y);
}

void FluidSolver::addSource(int x, int y)
{
    int index = getIndexOfCell(x, y);
    density[index] = source;
}

int FluidSolver::getIndexOfCell(int x, int y)
{
    ofVec2f cellPos = win2Cell(x, y);
    return idx(cellPos.x, cellPos.y);
}

//--------------------------------------------------------------
// Private
//--------------------------------------------------------------

void FluidSolver::addSource(std::vector<float>& x, const std::vector<float>& s)
{
    for (int i = 0; i < numCells; i++)
    {
        x[i] += s[i] * dt;
    }
}

void FluidSolver::setBound(int bound, std::vector<float>& x)
{
    int N = cellW;
    for(int i = 1; i <= N; i++)
    {
        x[idx(0  , i)] = bound == 1 ? -x[idx(1, i)] : x[idx(1,i)];
        x[idx(N+1, i)] = bound == 1 ? -x[idx(N, i)] : x[idx(N,i)];
        x[idx(i,   0)] = bound == 2 ? -x[idx(i, 1)] : x[idx(i,1)];
        x[idx(i, N+1)] = bound == 2 ? -x[idx(i, N)] : x[idx(i,N)];
    }

    x[idx(0,             0)] = 0.5 * (x[idx(1,           0)] + x[idx(0,           1)]);
    x[idx(0,       cellH+1)] = 0.5 * (x[idx(1,     cellH+1)] + x[idx(0,       cellH)]);
    x[idx(cellW+1,       0)] = 0.5 * (x[idx(cellW,       0)] + x[idx(cellW+1,     1)]);
    x[idx(cellW+1, cellH+1)] = 0.5 * (x[idx(cellW, cellH+1)] + x[idx(cellW+1, cellH)]);
}

void FluidSolver::linearSolve(int bound, std::vector<float>& x, std::vector<float>& x0, float a, float c)
{
    for (int k = 0; k < numItr; k++)
    {
        for(int i = 1; i <= cellW; i++)
        {
            for(int j = 1; j <= cellH; j++)
            {
                x[idx(i, j)] = (x0[idx(i, j)] + a * (x[idx(i-1, j)] + x[idx(i+1, j)] + x[idx(i, j-1)] + x[idx(i, j+1)])) / c;
            }
        }
        setBound(bound, x);
    }
}

void FluidSolver::diffuse(int bound, std::vector<float>& x, std::vector<float>& x0, float diff)
{
    float a = dt * diff * cellW * cellH;
    float c = 1 + 4 * a;
    linearSolve(bound, x, x0, a, c);
}

void FluidSolver::advect(int bound, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& u, std::vector<float>& v)
{
    float dt0X = dt * cellW;
    float dt0Y = dt * cellH;
    
    for(int i = 1; i <= cellW; i++)
    {
        for(int j = 1; j <= cellH; j++)
        {
            float x = i - dt0X * u[idx(i, j)];
            float y = j - dt0Y * v[idx(i, j)];
            
            if(x < 0.5)
                x = 0.5;
            if(x > cellW + 0.5)
                x = cellW + 0.5;
            
            int i0 = static_cast<int>(x);
            int i1 = i0 + 1;
            
            if(y < 0.5)
                y = 0.5;
            if(y > cellH + 0.5)
                y = cellH + 0.5;
            
            int j0 = static_cast<int>(y);
            int j1 = j0 + 1;
            
            float s1 = x - i0;
            float s0 = 1 - s1;
            float t1 = y - j0;
            float t0 = 1 - t1;
            
            d[idx(i, j)] = s0 * (t0 * d0[idx(i0, j0)] + t1 * d0[idx(i0, j1)]) +
                           s1 * (t0 * d0[idx(i1, j0)] + t1 * d0[idx(i1, j1)]);
        }
    }
    setBound(bound, d);
}

void FluidSolver::project(std::vector<float>& u, std::vector<float>& v, std::vector<float>& p, std::vector<float>& div)
{
    for (int i = 1; i <= cellW; i++)
    {
        for (int j = 1; j <= cellH; j++)
        {
            div[idx(i, j)] = -0.5 * (u[idx(i+1, j)] - u[idx(i-1, j)] + v[idx(i, j+1)] - v[idx(i, j-1)]) / cellW;
            p[idx(i, j)] = 0.0;
        }
    }
    setBound(0, div);
    setBound(0, p);
    
    linearSolve(0, p, div, 1, 4);
    
    for (int i = 1; i <= cellW; i++)
    {
        for (int j = 1; j <= cellH; j++)
        {
            u[idx(i, j)] -= 0.5 * cellW * (p[idx(i+1, j)] - p[idx(i-1, j)]);
            v[idx(i, j)] -= 0.5 * cellH * (p[idx(i, j+1)] - p[idx(i, j-1)]);
        }
    }
    setBound(1, u);
    setBound(2, v);
}

void FluidSolver::velocityStep()
{
    addSource(u, uPrev);
    addSource(v, vPrev);
    
    swap(uPrev, u);
    diffuse(1, u, uPrev, viscosity);
    swap(vPrev, v);
    diffuse(2, v, vPrev, viscosity);
    
    project(u, v, uPrev, vPrev);
    swap(uPrev, u);
    swap(vPrev, v);
    
    advect(1, u, uPrev, uPrev, vPrev);
    advect(2, v, vPrev, uPrev, vPrev);

    project(u, v, uPrev, vPrev);
}

void FluidSolver::densityStep()
{
    addSource(density, densityPrev);
    
    swap(densityPrev, density);
    diffuse(0, density, densityPrev, diff);
    
    swap(densityPrev, density);
    advect(0, density, densityPrev, u, v);
}

//--------------------------------------------------------------
// Utils
//--------------------------------------------------------------

inline void FluidSolver::swap(std::vector<float>& x0, std::vector<float>& x)
{
    x0.swap(x);
}

inline int FluidSolver::idx(int i, int j)
{
    return j * (cellW + 2) + i;
}

inline ofVec2f FluidSolver::win2Cell(int x, int y)
{
    float x_ = static_cast<float>(x);
    float y_ = static_cast<float>(y);
    
    float winW_ = static_cast<float>(winW);
    float winH_ = static_cast<float>(winH);
    
    float cellW_ = static_cast<float>(cellW);
    float cellH_ = static_cast<float>(cellH);
    
    float aspectX = cellW_ / winW_;
    float aspectY = cellH_ / winH_;
    
    float X = x_ * aspectX;
    float Y = y_ * aspectY;
    
    return ofVec2f(X, Y);
}
