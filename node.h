#pragma once

#include "bilinear_interpolation.h"

class Node
{
public:

    Node(const BilinearInterpolation& interpolation_);
    ~Node();

    bool isLeaf() const;

    void split(const BilinearInterpolation &bi1, const BilinearInterpolation &bi2);

    BilinearInterpolation interpolation;

    Node* pParent;
    Node* pChild1;
    Node* pChild2;

private:

    //disable copying
    Node(const Node&);
    Node& operator=(const Node&);
};
