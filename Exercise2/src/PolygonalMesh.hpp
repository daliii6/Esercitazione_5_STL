#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <list>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolygonalLibrary {

struct PolygonalMesh
{
    // Numero di celle di ogni tipo
    unsigned int NumCell0Ds;
    unsigned int NumCell1Ds;
    unsigned int NumCell2Ds;

    // Identificativi delle celle
    vector<unsigned int> Cell0DsId; ///< Cell0D id, size 1 x NumberCell0D
    vector<unsigned int> Cell1DsId; ///< Cell1D id, size 1 x NumberCell1D
    vector<unsigned int> Cell2DsId; ///< Cell2D id, size 1 x NumberCell2D

    // Coordinate dei punti (Cell0D)
    Eigen::MatrixXd Cell0DsCoordinates; // uso matrici per le funzioni di export di Vicini

    // Estremi dei segmenti (Cell1D): matrice 2xN con indici dei punti
    Eigen::MatrixXi Cell1DsExtrema;  ///< Cell1D vertices indices, size 2 x NumberCell1D (fromId,toId)

    // Vertici dei poligoni (Cell2D)
    vector<vector<unsigned int>> Cell2DsVertices; ///< Cell2D Vertices indices, size 1 x NumberCell2DVertices[NumberCell2D]

    // Lati dei poligoni (Cell2D)
    vector<vector<unsigned int>> Cell2DsEdges; ///< Cell2D Cell1D indices, size 1 x NumberCell2DEdges[NumberCell2D]

    // Marker associati alle celle: mappa da marker a lista di ID
    map<unsigned int, list<unsigned int>> Cell0DMarkers; ///< Cell0D markers
	map<unsigned int, list<unsigned int>> Cell1DMarkers; ///< Cell1D properties

};

}

