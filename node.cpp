#include "node.h"
#include <iostream>



Node::Node(const BilinearInterpolation &_interpolation) :
    interpolation(_interpolation),
    pParent(0),
    pChild1(0),
    pChild2(0)
{
}

Node::~Node()
{
    delete pChild1; pChild1 = NULL;
    delete pChild2; pChild2 = NULL;
}

bool Node::isLeaf() const
{
    return (NULL == pChild1) && (NULL == pChild2); 
}

void Node::split(const BilinearInterpolation &bi1,
                 const BilinearInterpolation &bi2)
{
    pChild1 = new Node(bi1);
    pChild2 = new Node(bi2);
    pChild1->pParent = this;
    pChild2->pParent = this;
}
