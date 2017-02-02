#ifndef JIAN_NUC2D_VIEW_VIEWSS
#define JIAN_NUC2D_VIEW_VIEWSS

#include<GL/glut.h>  

namespace jian {
namespace nuc2d {
namespace view {

inline void init() {
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glMatrixMode(GL_PROJECTION);
    glOrtho(-5, 5, -5, 5, 5, 15);
    glMatrixMode(GL_MODELVIEW);
    gluLookAt(0, 0, 10, 0, 0, 0, 0, 1, 0);
}

inline void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    glColor3f(1.0, 0, 0);
    glutWireTeapot(3);

    glFlush();
}

class ViewSS {
public:
    void operator ()(int argc, char* argv[]) {
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
        glutInitWindowPosition(0, 0);
        glutInitWindowSize(300, 300);

        glutCreateWindow("OpenGL 3D View");

        init();
        glutDisplayFunc(display);

        glutMainLoop();
    }
};


} // namespace view
} // namespace nuc2d
} // namespace jian

#endif

