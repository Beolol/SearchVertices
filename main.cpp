#include <QVector>
#include <QPair>
#include "FilterDate.h"
#define N 128
//#define SEGMENTS 4
using namespace std;

int main()
{
    size_t DegreePolinom = 6;

    qreal Y_to;

    QVector<qreal> func;
    for(size_t i = 0; i < N; i++)
    {
        Y_to=abs(qrand()/100);
        func.push_back( Y_to );
    }
    CFilterDate function(func, DegreePolinom, 15);
    function.GetAproxPoly();
    function.GetVertices(1, 128);

}
