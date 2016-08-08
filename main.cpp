#include <QVector>
#include <QPair>
#include <QCoreApplication>
#include "FilterDate.h"
#define N 128
#include <QDebug>


using namespace std;


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    size_t DegreePolinom = 6;

    qreal Y_to;

    QVector<qreal> func;
    for(size_t i = 0; i < N; i++)
    {
        Y_to=abs(qrand()/100);
        func.push_back( Y_to );
    }
    FilterData function(func, DegreePolinom, 15);
    auto poly = function.GetAproxPoly();
    auto vert = function.GetVertices(1, 128);
    for( int i = 0; i < poly.size(); i++ )
    {
        qDebug() << "x = "<< i << ", y = " + QString::number(poly[i]);
    }
    qDebug() << "Vertices:";
    for( auto v : vert )
    {
       qDebug() << "x = "<< v.first << " y = "<< v.second;
    }
    return a.exec();
}
