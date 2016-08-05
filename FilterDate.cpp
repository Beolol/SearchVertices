#include "FilterDate.h"
#include <QFile>
#include <QTextStream>



CFilterDate::CFilterDate(QVector<qreal> &date, const size_t Ndegree, const size_t SizeRWindow)
{
    size_t Ndeg = Ndegree;
    if(Ndeg>71) Ndeg = 70;
    int SizeRunningWindow = SizeRWindow;
    int size = date.size();
    int index;
    if((SizeRunningWindow%2)==0)
    {
        index = SizeRunningWindow/2;
        SizeRunningWindow++;
    }
    else
    {
        index = (SizeRunningWindow-1)/2;
    }
    for(int i = 0; i<index; i++)
    {
       QVector<qreal> sData;
       for(size_t j = 0; j< SizeRunningWindow; j++)
       {
           sData.push_back(date[j]);
       }
       PolySmoorthM(sData,Ndeg,i);
    }
    for(int i = index; i< (size - index ); i++)
    {
        QVector<qreal> sData;
        for(int j = (i-index); j< i+index+1; j++)
        {
            sData.push_back(date[j]);
        }
        PolySmoorthM(sData,Ndeg,index);
    }
    for(int i = size - index; i<size; i++)
    {
        QVector<qreal> sData;
        for(int j = size-SizeRunningWindow; j< size; j++)
        {
            sData.push_back(date[j]);
        }
        int point = index + (i-(size-index));
        PolySmoorthM(sData,Ndeg,point);
    }
    XTrAlgorithm(0,size-1);
}

void CFilterDate::PolySmoorthM(QVector<qreal> &sData,const size_t ndeg,const size_t index)
{
    QVector < QPair< qreal, qreal > > function;
    for( int i = 0; i<sData.size(); i++ )
    {
        function.push_back( QPair< qreal , qreal > ( i /*+ sData_length*/, sData[ i ] ) );
    }
    QScopedPointer<CPolinom> colPolinom(new CPolinom ( function.size(), ndeg + 1, function ));
    colPolinom->LeastSquareMethod();
    colPolinom->Gauss();
    value.push_back(colPolinom->ValueFunction(index));
    colPolinom->CalcDerivative();
    valueDev.push_back(colPolinom->ValueDerFunction(index));
    valueDev2.push_back(colPolinom->ValueSecondDerFunction(index));
    if(colPolinom->ValueSecondDerFunction(index) != 0)
    {
     qreal vCRadius = pow(1+pow(colPolinom->ValueDerFunction(index),2),3/2) / pow(pow(colPolinom->ValueSecondDerFunction(index),2),0.5);

     //qreal CRadius = (pow(1 + colPolinom->ValueDerFunction(index), 3/2)) / colPolinom->ValueSecondDerFunction(index);
     curveRadius.push_back( vCRadius );
    }
    else curveRadius.push_back( 0 );
}

QPair< QVector< qreal >, QVector< qreal > > CFilterDate::GetDebPoly
(
    QVector< qreal > &                      sData,
    const size_t                            ndeg,
    const size_t                            from
)
{
    QVector < QPair< qreal, qreal > > function;
    for( int i = 0; i<sData.size(); i++ )
    {
        function.push_back( QPair< qreal , qreal > ( i+from, sData[ i ] ) );
    }
    QScopedPointer< CPolinom > colPolinom( new CPolinom ( function.size(), ndeg + 1, function ) );
    colPolinom->LeastSquareMethod();
    colPolinom->Gauss();

    QVector< qreal > x;
    QVector< qreal > poly;

    for( int i = 0; i < sData.size(); ++i )
    {
        x.push_back( i + from );
        poly.push_back( colPolinom->ValueFunction( i ) );
    }

    return QPair< QVector< qreal >, QVector< qreal > >( x, poly );
}

void CFilterDate::Filter( QVector < qreal > & colDate, const size_t Ndegree ,const size_t SizeRunningWindow)
{
    size_t dSize = colDate.count();
    size_t from = 0;
    size_t to = qBound( from, from+SizeRunningWindow, dSize );
    while( from < dSize && to <= dSize )
    {
        ConstructionPolynomial(colDate,Ndegree,from,to);
        from += SizeRunningWindow;
        to = qBound( from, to + SizeRunningWindow, dSize );
    }
}

QVector<qreal> CFilterDate::GetAproxPoly()
{
//    QVector<qreal> VectorY;
//    for( int i =0; i < m_colValue.count(); i++ )
//    {
//        VectorY.push_back( m_colValue[i] );
//    }
//    return VectorY;

   return value;

}

void                                        CFilterDate::SavePolyToFile( const QString & filename )
{
    QFile workingFile( filename );
    if( workingFile.open( QFile::WriteOnly | QFile::Text ) )
    {
        QTextStream textBuffer( & workingFile );

        for( auto pair : aproxPoly )
        {
            textBuffer << pair.second << "\n";
        }

        workingFile.flush();
        workingFile.close();
    }
}

QVector< QPair< qreal, qreal > >              CFilterDate::GetVertices( qreal dInputA, qreal dInputB )
{
    QVector< QPair< qreal, qreal > > colMaxVertices;

    for( auto pair : vertices )
    {
        if( ( pair.first >= dInputA ) && ( pair.first <= dInputB ) )
        {
            colMaxVertices.push_back( pair );
#ifdef USE_STDOUT
            std::cout << pair.first << " : " << pair.second << std::endl;
#endif
        }
    }
    return colMaxVertices;
}

QVector< QPair< qreal, qreal > >              CFilterDate::GetMinVertices( qreal dInputA, qreal dInputB )
{
    QVector< QPair< qreal, qreal > > colMinVertices;

    for( auto pair : minVertices )
    {
        if( ( pair.first >= dInputA ) && ( pair.first <= dInputB ) )
        {
            colMinVertices.push_back( pair );
#ifdef USE_STDOUT
            std::cout << pair.first << " : " << pair.second << std::endl;
#endif
        }
    }
    return colMinVertices;
}

qreal CFilterDate::GetDerivative( const qreal &x )
{
    if( ( x >= 0) && (x < valueDev.count()))
        return valueDev[x];
//    if( ( x >= 0) && (x < m_colValueDerivative.count()))
//        return m_colValueDerivative[x];
    return 0;
}

qreal CFilterDate::GetSecondDrivative( const qreal &x )
{
    if( ( x >= 0) && (x < valueDev2.count()))
        return valueDev2[x];
//    if( ( x >= 0) && (x < m_colValueSecondDerivative.count()))
//        return m_colValueSecondDerivative[x];
    return 0;
}

qreal CFilterDate::GetCurveRadius( const qreal &x )
{
    if( ( x >= 0) && (x < curveRadius.count()))
        return curveRadius[x];
    return 0;
}

qreal CFilterDate::GetCurvative(const qreal &x)
{
    if( ( x >= 0) && (x < curvative.count()))
        return curvative[x];
    return 0;
}

void CFilterDate::XTrAlgorithm(const size_t &lBorder, const size_t &rBorder)
{
    size_t x_1 = lBorder;
    size_t x_2 = rBorder;
    for(size_t i = x_1; i < x_2; i++)
    {
        if(( valueDev[i] * valueDev[i+1] ) < 0 )
        {
          qreal Extrema_X =  i + ( valueDev[i]/(valueDev[i]-valueDev[i+1]) );
          qreal Extrema_Y = ( value[i+1]-value[i]) * Extrema_X + value[i]-i*(value[i+1]-value[i]);
          qreal D2V  = ( valueDev2[i+1]-valueDev2[i]) * Extrema_X + valueDev2[i]-i*(valueDev2[i+1]-valueDev2[i]);
          qreal Curvative = ( curveRadius[i+1] - curveRadius[i] ) * Extrema_X + curveRadius[i] - i * ( curveRadius[i+1] - curveRadius[i] );
          curvative.push_back(Curvative);
          if( D2V < 0 )
          {
              vertices.push_back( QPair< qreal, qreal >( Extrema_X, Extrema_Y ) );
          }
          if( D2V > 0)
          {
              minVertices.push_back( QPair< qreal, qreal >( Extrema_X, Extrema_Y ) );
          }
        }
    }

}

void CFilterDate::ConstructionPolynomial(QVector<qreal> &date, const size_t Ndegree, const size_t from, const size_t to)
{
    QVector < QPair< qreal, qreal > > function;
    size_t isize = abs((int)to-(int)from);
    for( size_t i = 0; (i + from <= to)&&(i+from<128); i++ )
    {
        function.push_back( QPair< qreal , qreal > ( i , date[ i + from ] ) );
    }
    QScopedPointer<CPolinom> colPolinom(new CPolinom ( isize, Ndegree + 1, function ));
    colPolinom->LeastSquareMethod();
    colPolinom->Gauss();
    for( size_t i = 0; i < isize; i++ )
    {
        aproxPoly.push_back( QPair< qreal , qreal >( i , colPolinom->ValueFunction( i ) ) );
    }
    colPolinom->CalcDerivative();
    colPolinom->SearchRoots( 0, isize - 1 );
    for( auto vert : colPolinom->ReturnMaxVertices() )
    {
        vertices.push_back( QPair< qreal , qreal > (vert.first + from, vert.second ));
    }
    for( auto vert : colPolinom->ReturnMinVertices() )
    {
        minVertices.push_back( QPair< qreal , qreal > (vert.first + from, vert.second ));
    }
    for( size_t i = 0 ; i < isize; i++)
    {
        valueDerivative.push_back( colPolinom->ValueDerFunction( i ) );
        valueSecondDerivative.push_back( colPolinom->ValueSecondDerFunction( i ) );
    }
}
