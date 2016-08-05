#include "Method.h"
#include <iostream>

#include <QDebug>

using namespace std;

CPolinom::CPolinom( size_t NSize, size_t Ndegree, const QVector< QPair< qreal, qreal > > & function )
{
    size = NSize;
    degree = Ndegree;

    vectorC.resize( degree );
    massiv.resize( degree );
    for( size_t count = 0; count < degree; count++ )
    {
        massiv[ count ].resize( degree + 1 );
    }
    for( auto pair : function )
    {
        inputX.push_back( pair.first );
        inputY.push_back( pair.second );
    }
}

qreal                                       CPolinom::ValueFunction( const qreal & x )
{
    qreal function = 0;
    for( size_t counter = 0; counter < degree; counter++ )
    {
        function += pow( x, counter ) * vectorC[ counter ];
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
    for( size_t i = 0; i < degree; i++ )
    {
        double a = massiv[ i ][ i ];
        for( size_t j = i + 1; j < degree; j++ )
        {
            double b = massiv[ j ][ i ];
            if( b != 0 )
            {
                for( size_t k = i; k < degree + 1; k++ )
                {
                    massiv[ j ][ k ] = massiv[ j ][ k ] * a / b
                            - massiv[ i ][ k ];
                }
            }

        }
    }
    //     PrintMassive();
    for( int i = degree - 1; i >= 0; i-- )
    {
        double summ = 0;
        for( size_t j = i + 1; j < degree; j++ )
        {
            summ += massiv[ i ][ j ] * vectorC[ j ];
        }
        summ = massiv[ i ][ degree ] - summ;
        if( massiv[ i ][ i ] == 0 )
        {
#ifdef USE_STDOUT
            std::cout << " Error " << std::endl;
#endif
            return;//@ToDo
        }
        vectorC[ i ] = summ / massiv[ i ][ i ];
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
    for( size_t counter = 1; counter < degree; counter++ )
    {
        derivative.push_back( vectorC[ counter ] * counter );
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
            if( ( degree - 1 ) == 1 )
            {
                newtonY = - derivative[ 0 ] / derivative[ 1 ];
            }
#ifdef USE_STDOUT
            std::cout <<  "Newton: " << newtonX << " y\"(x\') = " << newtonY << std::endl;
#endif
            if((newtonX>=x_1)&&(newtonX<=x_2))
            {
                if(( newtonY < 0))
                {
                    maxVertices.push_back( QPair< qreal, qreal >( newtonX, ValueFunction( newtonX ) ) );
                }
                else if( newtonY > 0)
                {
                    minVertices.push_back( QPair< qreal, qreal >( newtonX, ValueFunction( newtonX ) ) );
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

    if( ( degree - 2 ) == 1 )
    {
        auto valueSecondDerFunction = ValueSecondDerFunction( - derivative[ 0 ] / derivative[ 1 ] );
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
    for( size_t counter = 0; counter < degree; counter++ )
    {
        for( size_t j = 0; j < degree; j++ )
        {
            qreal sumA = 0;
            qreal sumB = 0;
            for( size_t k = 0; k < size; k++ )
            {
                sumA +=  pow( inputX[ k ], counter) * pow( inputX[ k ], j); //pow( m_colInputX[k], counter + j );//
                sumB +=  inputY[ k ] * pow( inputX[ k ], counter );
            }
            massiv[ counter ][ j ] = sumA;
            massiv[ counter ][ degree ] = sumB;
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
    for( size_t counter = 0; counter < degree - 1; counter++ )
    {
        function += pow( x, counter ) * derivative[ counter ];
    }
    return function;
}

qreal                                       CPolinom::ValueSecondDerFunction( const qreal & x )
{
    qreal function = 0;
    for( size_t counter = 1; counter < degree - 1; counter++ )
    {
        function += pow( x, counter - 1 ) * derivative[ counter ] * counter;
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
    for( size_t counter = 1; counter < degree - 1; counter++ )
    {
        f_derivative += pow( x1, counter - 1 ) * derivative[ counter ] * counter;
    }
    if( f_derivative != 0 )
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
    for( size_t counter = 0; counter < size; counter++ )
    {
        result[ counter ].first = counter;
        result[ counter ].second = ValueDerFunction( counter );
    }
    return result;
}

QVector< QPair< qreal, qreal > >            CPolinom::ReturnMaxVertices()
{
    return maxVertices;
}

QVector< QPair< qreal, qreal > >            CPolinom::ReturnMinVertices()
{
    return minVertices;
}

void                                        CPolinom::Frebenius()
{
    QVector< QVector< qreal > > matrix( degree - 1 );
    for( size_t count = 0; count < degree - 1; count++ )
    {
        matrix[ count ].resize( degree - 1 );
    }
    for( size_t counter = 0; counter < degree - 1; counter++ )
    {
        for( size_t j = 0; j < degree - 1; j++ )
        {
            if( counter == 0 )
            {
                qreal xx = vectorC[ j ];
                qreal yy = vectorC[ degree - 1 ];
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
