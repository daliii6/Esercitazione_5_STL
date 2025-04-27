#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolygonalLibrary
{

bool ImportMesh(PolygonalMesh& mesh)
{
    if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;
    
    return true;
}

// ***************************************************************************
bool ImportCell0Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell0Ds.csv");
    if(file.fail()) return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();
    listLines.pop_front();  // Rimuove l'header

    mesh.NumCell0Ds = listLines.size();
    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (string& line : listLines)
    {
        replace(line.begin(), line.end(), ';', ' '); // Cambia i separatori, per poter usare quanto fatto in Exercise1
        istringstream converter(line);

        unsigned int id, marker;
        Vector2d coord;

        converter >> id >> marker >> mesh.Cell0DsCoordinates(0, id) >> mesh.Cell0DsCoordinates(1, id);
        /*Nella prima riga di Cell0DsCoordinates, abbiamo la coordinata x del punto.
        Nella seconda riga la coordinata y. L'id corrisponde alla colonna.*/
        
        mesh.Cell0DsId.push_back(id);

        // Aggiungo i marker nella mappa
        if(marker != 0) //solo quelli laterali
		{
			
			const auto it = mesh.Cell0DMarkers.find(marker);
			if(it == mesh.Cell0DMarkers.end())  // se il marker esiste già
                mesh.Cell0DMarkers.insert({marker, {id}}); // aggiungo l'id alla lista
			else
                //mesh.Cell0DMarkers[marker].push_back(id);
                it->second.push_back(id);
            /*Se il marker non è stato ancora aggiunto, 
            si inserisce una nuova voce nella mappa con il marker 
            come chiave e una lista contenente solo id come valore*/
			
		}
    }

    return true;
}



// ***************************************************************************
bool ImportCell1Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell1Ds.csv");
    if(file.fail()) return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();
    listLines.pop_front();  // Rimuove l'header

    mesh.NumCell1Ds = listLines.size();
    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

    for (string& line : listLines)
    {
        replace(line.begin(), line.end(), ';', ' ');
        istringstream converter(line);

        unsigned int id, marker;
        Vector2i vertices;

        converter >> id >> marker >> mesh.Cell1DsExtrema(0, id) >> mesh.Cell1DsExtrema(1, id);
        mesh.Cell1DsId.push_back(id);

        /// Memorizza i marker
        if(marker != 0)
        {
            const auto it = mesh.Cell1DMarkers.find(marker);
            if(it == mesh.Cell1DMarkers.end())
            {
                mesh.Cell1DMarkers.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell1Ds[marker].push_back(id);
                it->second.push_back(id);
            }
        }
    }

    return true;
}


// ***************************************************************************
bool ImportCell2Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell2Ds.csv");
    if (file.fail()) return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();
    listLines.pop_front();  // Rimuove l'header

    mesh.NumCell2Ds = listLines.size();
    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (string& line : listLines)
    {
        replace(line.begin(), line.end(), ';', ' '); // Cambia i separatori
        istringstream converter(line);

        unsigned int id, marker, numVertices, numEdges;

        // Leggo id, marker, numero di vertici
        converter >> id >> marker >> numVertices;

        // Leggo i vertici e li inserisco in un vettore 
        vector<unsigned int> vertices(numVertices);  
        for (unsigned int i = 0; i < numVertices; i++)
            converter >> vertices[i];

        // Leggo numero di spigoli
        converter >> numEdges;

        // Leggo gli spigoli
        vector<unsigned int> edges(numEdges);
        for (unsigned int i = 0; i < numEdges; i++)
            converter >> edges[i];

        // Memorizza l’id
        mesh.Cell2DsId.push_back(id);

        // Inserisci nei vettori dei vertici e degli spigoli
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsEdges.push_back(edges);
    
    }

    return true;
}


}