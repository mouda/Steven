/*
  File: ULAGENT.h
  Brief:1.each object,ULAGENT is one machine/node.
        2.ULAGENT only called/managed by ULSA3b

  Author: Steven
  Date: 2012/08/07
  Latest Progress:

  Abbreviation:

*/


#include<vector>

#include "ULAGENT.h"
using namespace std;

ULAGENT::ULAGENT()
{existed = false;}

ULAGENT::ULAGENT(const ULAGENT & ref)
{
    existed = ref.existed;
    nodeIndex = ref.nodeIndex;
    locX = ref.locX;
    locY = ref.locY;
    ptrHead= ref.ptrHead;
    power = 1e-11;//set as a value which is close to zero
    T2ndtier=0;
}

/*
 Purpose: construct a background info for a single machine/agent/node
*/
void ULAGENT::aryConstructor(int inputIndex, float X, float Y)
{
    existed = true;
    nodeIndex = inputIndex;// please Note nodeIndex := nodeName
    locX = X;
    locY = Y;
    ptrHead = NULL;
    dibSQ=locX*locX+locY+locY*locY;
}
