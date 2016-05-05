#include "Method.h"
#include <iostream>

#include <QDebug>

using namespace std;

CPolinom::CPolinom( size_t NSize, size_t Ndegree, const QVector< QPair< qreal, qreal > > & function )
{
    m_iSize = NSize;
    m_iDegree = Ndegree;

    m_colVectorC.resize( m_iDegree );
    m_colMassiv.resize( m_iDegree );
    for( size_t count = 0; count < m_iDegree; count++ )
    {
        m_colMassiv[ count ].resize( m_iDegree + 1 );
    }
    for( auto pair : function )
    {
        m_colInputX.push_back( pair.first );
        m_colInputY.push_back( pair.second );
    }
}

qreal                                       CPolinom::ValueFunction( const qreal & x )
{
    qreal function = 0;
    for( size_t counter = 0; counter < m_iDegree; counter++ )
    {
        function += pow( x, counter ) * m_colVectorC[ counter ];
    }
    return function;
}

void                                        CPolinom::PrintMassive()
{
#ifdef USE_STDOUT
    std::cout << "matrix " << m_iDegree << "x" << m_iDegree << ":" << std::endl;
    for( size_t i = 0; i < m_iDegree + 1; i++ )
    {
        for( size_t j = 0; j < m_iDegree; j++ )
        {
            std::cout << m_colMassiv[ j ][ i ] << " ";
        }
        std::cout << std::endl;
    }
#endif
}


void                                        CPolinom::Gauss()
{
#ifdef USE_STDOUT
    std::cout<<"input Gauss"<<std::endl;
#endif
    //     PrintMassive();
    for( size_t i = 0; i < m_iDegree; i++ )
    {
        double a = m_colMassiv[ i ][ i ];
        for( size_t j = i + 1; j < m_iDegree; j++ )
        {
            double b = m_colMassiv[ j ][ i ];
            if( b != 0 )
            {
                for( size_t k = i; k < m_iDegree + 1; k++ )
                {
                    m_colMassiv[ j ][ k ] = m_colMassiv[ j ][ k ] * a / b
                            - m_colMassiv[ i ][ k ];
                }
            }

        }
    }
    //     PrintMassive();
    for( int i = m_iDegree - 1; i >= 0; i-- )
    {
        double summ = 0;
        for( size_t j = i + 1; j < m_iDegree; j++ )
        {
            summ += m_colMassiv[ i ][ j ] * m_colVectorC[ j ];
        }
        summ = m_colMassiv[ i ][ m_iDegree ] - summ;
        if( m_colMassiv[ i ][ i ] == 0 )
        {
#ifdef USE_STDOUT
            std::cout << " Error " << std::endl;
#endif
            return;//@ToDo
        }
        m_colVectorC[ i ] = summ / m_colMassiv[ i ][ i ];
    }

#ifdef USE_STDOUT
    std::cout << " y = ";
    for( size_t i = 0; i < m_iDegree; i++ )
    {
        if( ( i != 0) && ( m_colVectorC[ i ] > 0 ) )
        {
            std::cout << "+";
        }
        std::cout << m_colVectorC[ i ] << "*x^" << i << " ";
    }
    std::cout <<"out Gauss"<<std::endl;
#endif
}

void                                        CPolinom::CalcDerivative()
{
#ifdef USE_STDOUT
    std::cout <<"input Derivative"<<std::endl;
#endif
    for( size_t counter = 1; counter < m_iDegree; counter++ )
    {
        m_colDerivative.push_back( m_colVectorC[ counter ] * counter );
    }
#ifdef USE_STDOUT
    std::cout << " y'= ";
    for( size_t counter = 0; counter < m_iDegree - 1; counter++ )
    {
        if( ( counter !=0 ) && ( m_colDerivative[ counter ] > 0 ) )
        {
            std::cout << "+";
        }
        std::cout << m_colDerivative[ counter ] << "*x^" << counter << " ";
    }
    std::cout << std::endl;
    std::cout <<"output Derivative"<<std::endl;
#endif
}


void                                        CPolinom::SearchRoots( const qreal & lBorder, const qreal & section )
{
#ifdef USE_STDOUT
    std::cout <<"input SearchRoots"<<std::endl;
#endif
    qreal eror = 0.0001;
    qreal x_1 = lBorder;
    qreal div = 128 - section;
    if( div == 0 )
    {
        div = 1;
    }
    qreal x_i = 1.0 / div;
    qreal x_2 = lBorder + x_i;
    bool flag = false;
    qreal counter = 0;
    while( ( counter <= section) )
    {
        if( ValueDerFunction( x_1 ) * ValueDerFunction( x_2 ) < 0 )
        {
            if( ValueSecondDerFunction( x_1 ) * ValueSecondDerFunction( x_2 ) > 0 )
            {
                flag = true;
            }
            //            else x_1 -= x_i;
            qreal newtonX1 = x_1;
            qreal newtonX2 = x_2;
            auto methodSecants = MethodSecants( newtonX1, newtonX2, eror);
#ifdef USE_STDOUT
            std::cout << "MethodSecants: " << methodSecants << std::endl;
#endif
        }
        if( flag )
        {
            qreal newtonX1 = x_1;
            qreal newtonX2 = x_2;
#ifdef USE_STDOUT
            std::cout << newtonX1<< " "<< newtonX2;
#endif
            qreal newtonX = Newton( newtonX2, newtonX1, eror );
            qreal newtonY = ValueSecondDerFunction( newtonX );
            if( ( m_iDegree - 1 ) == 1 )
            {
                newtonY = - m_colDerivative[ 0 ] / m_colDerivative[ 1 ];
            }
#ifdef USE_STDOUT
            std::cout <<  "Newton: " << newtonX << " y\"(x\') = " << newtonY << std::endl;
#endif
            if((newtonX>=x_1)&&(newtonX<=x_2))
            {
                if(( newtonY < 0))
                {
                    m_colMaxVertices.push_back( QPair< qreal, qreal >( newtonX, ValueFunction( newtonX ) ) );
                }
                else if( newtonY > 0)
                {
                    m_colMinVertices.push_back( QPair< qreal, qreal >( newtonX, ValueFunction( newtonX ) ) );
                }
            }
            flag = false;
            x_1 = x_2;
            x_2 +=x_i;
        }
        else
        {
            x_1 += x_i;
            x_2 += x_i;
        }
        counter+=x_i;
    }

    if( ( m_iDegree - 2 ) == 1 )
    {
        auto valueSecondDerFunction = ValueSecondDerFunction( - m_colDerivative[ 0 ] / m_colDerivative[ 1 ] );
#ifdef USE_STDOUT
        std::cout << - m_colDerivative[ 0 ] / m_colDerivative[ 1 ]
                << " y\"(x\') = "
                <<  valueSecondDerFunction << std::endl;
#endif
    }
#ifdef USE_STDOUT
    std::cout <<"output SearchRoots"<<std::endl;
#endif
}


void                                        CPolinom::LeastSquareMethod()
{
#ifdef USE_STDOUT
    std::cout <<"input LeastSquareMethod"<<std::endl;
#endif
    for( size_t counter = 0; counter < m_iDegree; counter++ )
    {
        for( size_t j = 0; j < m_iDegree; j++ )
        {
            qreal sumA = 0;
            qreal sumB = 0;
            for( size_t k = 0; k < m_iSize; k++ )
            {
                sumA +=  pow( m_colInputX[ k ], counter) * pow( m_colInputX[ k ], j); //pow( m_colInputX[k], counter + j );//
                sumB +=  m_colInputY[ k ] * pow( m_colInputX[ k ], counter );
            }
            m_colMassiv[ counter ][ j ] = sumA;
            m_colMassiv[ counter ][ m_iDegree ] = sumB;
        }
    }
#ifdef USE_STDOUT
    std::cout <<"output LeastSquareMethod"<<std::endl;
#endif
}


qreal                                       CPolinom::MethodSecants( qreal & x1, qreal & x2, const qreal & e )
{
#ifdef USE_STDOUT
    std::cout <<"input MethodSecants"<<std::endl;
#endif
    if( fabs( x2 - x1 ) < e)
    {
        return x2;
    }
    x2 = ( x1
           -
           ( ( ValueDerFunction( x1 ) ) * ( x2 - x1 ) )
           /
           ( ValueDerFunction( x2 ) - ValueDerFunction( x1 ) )
           );
    return MethodSecants( x2, x1, e );

#ifdef USE_STDOUT
    std::cout <<"output MethodSecants"<<std::endl;
#endif
}
qreal                                       CPolinom::ValueDerFunction( const qreal & x )
{
    qreal function = 0;
    for( size_t counter = 0; counter < m_iDegree - 1; counter++ )
    {
        function += pow( x, counter ) * m_colDerivative[ counter ];
    }
    return function;
}

qreal                                       CPolinom::ValueSecondDerFunction( const qreal & x )
{
    qreal function = 0;
    for( size_t counter = 1; counter < m_iDegree - 1; counter++ )
    {
        function += pow( x, counter - 1 ) * m_colDerivative[ counter ] * counter;
    }
    return function;
}

qreal                                       CPolinom::Newton( qreal & x1, qreal & x2, const qreal & e )
{
#ifdef USE_STDOUT
    std::cout <<"input Newton"<<std::endl;
#endif
    if( fabs( x2 - x1 ) < e )
    {
#ifdef USE_STDOUT
        std::cout <<"output Newton1"<<std::endl;
#endif
        return x2;
    }
    qreal f_derivative = 0;
    for( size_t counter = 1; counter < m_iDegree - 1; counter++ )
    {
        f_derivative += pow( x1, counter - 1 ) * m_colDerivative[ counter ] * counter;
    }
    if( f_derivative )
    {
        x2 = x1;
        x1 = x1 - ValueDerFunction( x1 ) / f_derivative;
        return Newton( x1, x2, e );
    }
#ifdef USE_STDOUT
    std::cout <<"output Newton2"<<std::endl;
#endif
    return 0;
}

QVector< QPair< qreal, qreal > >            CPolinom::ReturnFunction()
{
    QVector < QPair< qreal, qreal > > result;
    for( size_t counter = 0; counter < m_iSize; counter++ )
    {
        result[ counter ].first = counter;
        result[ counter ].second = ValueDerFunction( counter );
    }
    return result;
}

QVector< QPair< qreal, qreal > >            CPolinom::ReturnMaxVertices()
{
    return m_colMaxVertices;
}

QVector< QPair< qreal, qreal > >            CPolinom::ReturnMinVertices()
{
    return m_colMinVertices;
}

void                                        CPolinom::Frebenius()
{
    QVector< QVector< qreal > > matrix( m_iDegree - 1 );
    for( size_t count = 0; count < m_iDegree - 1; count++ )
    {
        matrix[ count ].resize( m_iDegree - 1 );
    }
    for( size_t counter = 0; counter < m_iDegree - 1; counter++ )
    {
        for( size_t j = 0; j < m_iDegree - 1; j++ )
        {
            if( counter == 0 )
            {
                qreal xx = m_colVectorC[ j ];
                qreal yy = m_colVectorC[ m_iDegree - 1 ];
                matrix[ counter ][ j ] = - xx / yy;
            }
            else if( ( counter - j == 1 ) )
            {
                matrix[ counter ][ j ] = -1;
            }
            else
            {
                matrix[ counter ][ j ] = 0;
            }
#ifdef USE_STDOUT
            std::cout << " " << matrix[ counter ][ j ];
#endif
        }
#ifdef USE_STDOUT
        std::cout << std::endl;
#endif
    }
}
