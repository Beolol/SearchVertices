#include "FilterDate.h"
#include <iostream>
#include <QFile>
#include <QTextStream>



CFilterDate::CFilterDate(QVector<qreal> &date, const size_t Ndegree, const size_t SizeRWindow)
{
    size_t Ndeg = Ndegree;
    if(Ndeg>71) Ndeg = 70;
    size_t SizeRunningWindow = SizeRWindow;
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

void                                        CFilterDate::PolySmoorthM(QVector<qreal> &sData,const size_t ndeg,const size_t index)
{
    QVector < QPair< qreal, qreal > > function;
    for( int i = 0; i<sData.size(); i++ )
    {
        function.push_back( QPair< qreal , qreal > ( i /*+ sData_length*/, sData[ i ] ) );
    }
    QScopedPointer<CPolinom> colPolinom(new CPolinom ( function.size(), ndeg + 1, function ));
    colPolinom->LeastSquareMethod();
    colPolinom->Gauss();
    m_colValue.push_back(colPolinom->ValueFunction(index));
    colPolinom->CalcDerivative();
    m_colValueDev.push_back(colPolinom->ValueDerFunction(index));
    m_colValueDev2.push_back(colPolinom->ValueSecondDerFunction(index));
    if(colPolinom->ValueSecondDerFunction(index) != 0)
    {
     qreal vCRadius = pow(1+pow(colPolinom->ValueDerFunction(index),2),3/2) / pow(pow(colPolinom->ValueSecondDerFunction(index),2),0.5);

     //qreal CRadius = (pow(1 + colPolinom->ValueDerFunction(index), 3/2)) / colPolinom->ValueSecondDerFunction(index);
     m_colCurveRadius.push_back( vCRadius );
    }
    else m_colCurveRadius.push_back( 0 );
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

void                                        CFilterDate::Filter( QVector < qreal > & colDate, const size_t Ndegree ,const size_t SizeRunningWindow)
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

   return m_colValue;

}

void                                        CFilterDate::SavePolyToFile( const QString & filename )
{
    QFile workingFile( filename );
    if( workingFile.open( QFile::WriteOnly | QFile::Text ) )
    {
        QTextStream textBuffer( & workingFile );

        for( auto pair : m_colAproxPoly )
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

    for( auto pair : m_colVertices )
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

    for( auto pair : m_colMinVertices )
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
    if( ( x >= 0) && (x < m_colValueDev.count()))
        return m_colValueDev[x];
//    if( ( x >= 0) && (x < m_colValueDerivative.count()))
//        return m_colValueDerivative[x];
    return 0;
}

qreal CFilterDate::GetSecondDrivative( const qreal &x )
{
    if( ( x >= 0) && (x < m_colValueDev2.count()))
        return m_colValueDev2[x];
//    if( ( x >= 0) && (x < m_colValueSecondDerivative.count()))
//        return m_colValueSecondDerivative[x];
    return 0;
}

qreal CFilterDate::GetCurveRadius( const qreal &x )
{
    if( ( x >= 0) && (x < m_colCurveRadius.count()))
        return m_colCurveRadius[x];
    return 0;
}

qreal CFilterDate::GetCurvative(const qreal &x)
{
    if( ( x >= 0) && (x < m_colCurvative.count()))
        return m_colCurvative[x];
    return 0;
}

void CFilterDate::XTrAlgorithm(const size_t &lBorder, const size_t &rBorder)
{
    size_t x_1 = lBorder;
    size_t x_2 = rBorder;
    for(size_t i = x_1; i < x_2; i++)
    {
        if(( m_colValueDev[i] * m_colValueDev[i+1] ) < 0 )
        {
          qreal Extrema_X =  i + ( m_colValueDev[i]/(m_colValueDev[i]-m_colValueDev[i+1]) );
          qreal Extrema_Y = ( m_colValue[i+1]-m_colValue[i]) * Extrema_X + m_colValue[i]-i*(m_colValue[i+1]-m_colValue[i]);
          qreal D2V  = ( m_colValueDev2[i+1]-m_colValueDev2[i]) * Extrema_X + m_colValueDev2[i]-i*(m_colValueDev2[i+1]-m_colValueDev2[i]);
          qreal Curvative = ( m_colCurveRadius[i+1] - m_colCurveRadius[i] ) * Extrema_X + m_colCurveRadius[i] - i * ( m_colCurveRadius[i+1] - m_colCurveRadius[i] );
          m_colCurvative.push_back(Curvative);
          if( D2V < 0 )
          {
              m_colVertices.push_back( QPair< qreal, qreal >( Extrema_X, Extrema_Y ) );
          }
          if( D2V > 0)
          {
              m_colMinVertices.push_back( QPair< qreal, qreal >( Extrema_X, Extrema_Y ) );
          }
        }
    }

}

void CFilterDate::ConstructionPolynomial(QVector<qreal> &date, const size_t Ndegree, const size_t from, const size_t to)
{
    QVector < QPair< qreal, qreal > > function;
    size_t isize = abs(to-from);
    for( size_t i = 0; (i + from <= to)&&(i+from<128); i++ )
    {
        function.push_back( QPair< qreal , qreal > ( i , date[ i + from ] ) );
    }
    QScopedPointer<CPolinom> colPolinom(new CPolinom ( isize, Ndegree + 1, function ));
    colPolinom->LeastSquareMethod();
    colPolinom->Gauss();
    for( size_t i = 0; i < isize; i++ )
    {
        m_colAproxPoly.push_back( QPair< qreal , qreal >( i , colPolinom->ValueFunction( i ) ) );
    }
    colPolinom->CalcDerivative();
    colPolinom->SearchRoots( 0, isize - 1 );
    for( auto vert : colPolinom->ReturnMaxVertices() )
    {
        m_colVertices.push_back( QPair< qreal , qreal > (vert.first + from, vert.second ));
    }
    for( auto vert : colPolinom->ReturnMinVertices() )
    {
        m_colMinVertices.push_back( QPair< qreal , qreal > (vert.first + from, vert.second ));
    }
    for( size_t i = 0 ; i < isize; i++)
    {
        m_colValueDerivative.push_back( colPolinom->ValueDerFunction( i ) );
        m_colValueSecondDerivative.push_back( colPolinom->ValueSecondDerFunction( i ) );
    }
}
