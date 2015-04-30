//
// Created by Guillaume Rebmann on 28/04/15.
//

#include "BinomialConvolver.h"
/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file BinomialConvolver.ih
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2011/07/06
 *
 * Implementation of inline methods defined in BinomialConvolver.h
 *
 * This file is part of the DGtal library.
 */


//////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
//////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION of inline methods.
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// ----------------------- Standard services ------------------------------

//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>::~BinomialConvolver()
{
}
//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::BinomialConvolver( unsigned int n )
{
    setSize( n ); //std::cout << "Hello the world" << std::endl;
}
//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
void
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::setSize( unsigned int n )
{
    myN = n;
}
//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
unsigned int
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::size() const
{
    return myN;
}
//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
int
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::index( const ConstIteratorOnPoints& it ) const
{
    typename std::map<ConstIteratorOnPoints,int>::const_iterator
            map_it = myMapIt2Idx.find( it );
    if ( map_it != myMapIt2Idx.end() )
        return map_it->second;
    ASSERT( false );
    return 0;
}

//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
void
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::init( const double h,
        const ConstIteratorOnPoints& itb,
        const ConstIteratorOnPoints& ite,
        const bool isClosed )
{
    myMapIt2Idx.clear();
    myH = h;
    myBegin = itb;
    myEnd = ite;
    unsigned int aSize = 0;
    for ( ConstIteratorOnPoints it = itb; it != ite; ++it )
    {
        myMapIt2Idx[ it ] = aSize;
        ++aSize;
    }
    myX.init( aSize, 0, isClosed, 0.0 );
    myY.init( aSize, 0, isClosed, 0.0 );
    aSize = 0;
    for ( ConstIteratorOnPoints it = itb; it != ite; ++it, ++aSize )
    {
/*      myX[ size ] = it->operator[]( 0 );
      myY[ size ] = it->operator[]( 1 );*/
// TRIS ConstIterator may have no -> operator
        Point p(*it);
        myX[ aSize ] = p[0];
        myY[ aSize ] = p[1];
    }
    Signal<double> G = Signal<double>::G2n( myN );
    myX = myX * G;
    myY = myY * G;
    myDX = myX * Signal<double>::Delta();
    myDY = myY * Signal<double>::Delta();
    myDDX = myDX * Signal<double>::Delta();
    myDDY = myDY * Signal<double>::Delta();
}

//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
std::pair<TValue,TValue>
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::x( int i ) const
{
    return std::make_pair( myX[ i ], myY[ i ] );
}
//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
std::pair<TValue,TValue>
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::dx( int i ) const
{
    return std::make_pair( myDX[ i ], myDY[ i ] );
}
//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
std::pair<TValue,TValue>
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::d2x( int i ) const
{
    return std::make_pair( myDDX[ i ], myDDY[ i ] );
}
//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
std::pair<TValue,TValue>
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::tangent( int i ) const
{
    Value n = sqrt( myDX[ i ] * myDX[ i ] +
                    myDY[ i ] * myDY[ i ] );
    return std::make_pair( -myDX[ i ] / n, -myDY[ i ] / n );
}
//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
TValue
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::curvature( int i ) const
{
    Value denom = pow( myDX[ i ] * myDX[ i ] + myDY[ i ] * myDY[ i ], 1.5 );
    return ( denom != TValue( 0.0 ) )
           ? ( myDDX[ i ] * myDY[ i ] - myDDY[ i ] * myDX[ i ] ) / denom / myH
           : TValue( 0.0 );
}

/**
   @return the suggested size for the binomial convolver as
    ceil( d / pow( h, 1.0/3.0 ) ), with d the diameter of the
   contour.
*/
//-----------------------------------------------------------------------------
template <typename TConstIteratorOnPoints, typename TValue>
inline
unsigned int
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>
::suggestedSize( const double h,
                 const ConstIteratorOnPoints& itb,
                 const ConstIteratorOnPoints& ite )
{
    Point p(*itb);
    TValue xmin = p[ 0 ];
    TValue ymin = p[ 1 ];
    TValue xmax = p[ 0 ];
    TValue ymax = p[ 1 ];
    for ( ConstIteratorOnPoints it = itb; it != ite; ++it )
    {
/*      TValue x = it->operator[]( 0 );
      TValue y = it->operator[]( 1 );*/
// TRIS ConstIterator may have no -> operator
        Point pp(*it);
        TValue x = pp[0];
        TValue y = pp[1];
        if ( x < xmin ) xmin = x;
        if ( x > xmax ) xmax = x;
        if ( y < ymin ) ymin = y;
        if ( y > ymax ) ymax = y;
    }
    TValue diameter = ( xmax - xmin ) > ( ymax - ymin )
                      ? ( xmax - xmin )
                      : ( ymax - ymin );
//  return (unsigned int) ceil( 0.5 / pow( h / diameter, 4.0/3.0 ) );
//TRIS  (diameter*h is the diameter of the shape)
    return (unsigned int) ceil( diameter / pow( h, 1.0/3.0 ) );
}

///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param out the output stream where the object is written.
 */
template <typename TConstIteratorOnPoints, typename TValue>
inline
void
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>::selfDisplay ( std::ostream & out ) const
{
    out << "[BinomialConvolver]";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
double
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::getMyH() {
    return myBC.myH;
}

template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
unsigned int
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::getMyN() {
    return myBC.myN;
}

template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
DGtal::Signal<double>
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::getMyX() {
    return myBC.myX;
}

template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
DGtal::Signal<double>
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::getMyY() {
    return myBC.myY;
}

template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
DGtal::Signal<double>
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::getMyDX() {
    return myBC.myDX;
}
template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
DGtal::Signal<double>
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::getMyDY() {
    return myBC.myDY;
}

template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
DGtal::Signal<double>
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::getMyDDX() {
    return myBC.myDDX;
}
template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
DGtal::Signal<double>
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::getMyDDY() {
    return myBC.myDDY;
}


/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
template <typename TConstIteratorOnPoints, typename TValue>
inline
bool
DGtal::BinomialConvolver<TConstIteratorOnPoints,TValue>::isValid() const
{
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// TangentFromBinomialConvolverFunctor<,TBinomialConvolver,TRealPoint>
//-----------------------------------------------------------------------------
template <typename TBinomialConvolver, typename TRealPoint>
inline
typename DGtal::TangentFromBinomialConvolverFunctor<TBinomialConvolver,TRealPoint>::Value
DGtal::TangentFromBinomialConvolverFunctor<TBinomialConvolver,TRealPoint>
::operator()( const BinomialConvolver & bc,
              const ConstIteratorOnPoints & it ) const
{
    int index = bc.index( it );
    std::pair<SignalValue,SignalValue> v = bc.tangent( index );
    return RealPoint( v.first, v.second );
}

///////////////////////////////////////////////////////////////////////////////
// CurvatureFromBinomialConvolverFunctor<,TBinomialConvolver,TRealPoint>
//-----------------------------------------------------------------------------
template <typename TBinomialConvolver, typename TReal>
inline
typename DGtal::CurvatureFromBinomialConvolverFunctor<TBinomialConvolver,TReal>::Value
DGtal::CurvatureFromBinomialConvolverFunctor<TBinomialConvolver,TReal>
::operator()( const BinomialConvolver & bc,
              const ConstIteratorOnPoints & it ) const
{
    int index = bc.index( it );
    Value v = bc.curvature( index );
    return v;
}

///////////////////////////////////////////////////////////////////////////////
// class BinomialConvolverEstimator <TBinomialConvolver,TBinomialConvolverFunctor>
//-----------------------------------------------------------------------------
template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::BinomialConvolverEstimator( unsigned int n,
                              const BinomialConvolverFunctor & f )
        : myBC( n ), myFunctor( f )
{
}
//-----------------------------------------------------------------------------
template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
void
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::init( const double h,
        const ConstIterator & itb,
        const ConstIterator & ite,
        const bool isClosed )
{
    if ( myBC.size() == 0 )
        myBC.setSize( myBC.suggestedSize( h, itb, ite ) );
    myBC.init( h, itb, ite, isClosed );
}



//-----------------------------------------------------------------------------
template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
inline
typename DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>::Quantity
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::eval( const ConstIterator& it )
{
    return myFunctor( myBC, it );
}
//-----------------------------------------------------------------------------
template <typename TBinomialConvolver, typename TBinomialConvolverFunctor>
template <typename OutputIterator>
inline
OutputIterator
DGtal::BinomialConvolverEstimator<TBinomialConvolver,TBinomialConvolverFunctor>
::eval( const ConstIterator& itb,
        const ConstIterator& ite,
        OutputIterator result )
{
    for ( std::vector<Z2i::Point>::const_iterator it = itb; it != ite; ++it ){
        int i = 0;
        typename std::map<ConstIteratorOnPoints,int>::const_iterator
                map_it = myBC.myMapIt2Idx.find( it );
        if ( map_it != myBC.myMapIt2Idx.end() )
            i = map_it->second;
        double denom = pow( myBC.myDX[ i ] * myBC.myDX[ i ] + myBC.myDY[ i ] * myBC.myDY[ i ], 1.5 );
        double v = ( denom != double( 0.0 ) )? ( myBC.myDDX[ i ] * myBC.myDY[ i ] - myBC.myDDY[ i ] * myBC.myDX[ i ] ) / denom / myBC.myH : double( 0.0 );
        //std::cout << "Value: " << v << std::endl;
        *result++ = v;
    }
    return result;
}


///////////////////////////////////////////////////////////////////////////////
// Implementation of inline functions                                        //

template <typename TConstIteratorOnPoints, typename TValue>
inline
std::ostream&
DGtal::operator<<
        ( std::ostream & out,
          const BinomialConvolver<TConstIteratorOnPoints,TValue> & object )
{
    object.selfDisplay( out );
    return out;
}

//                                                                           //
///////////////////////////////////////////////////////////////////////////////


