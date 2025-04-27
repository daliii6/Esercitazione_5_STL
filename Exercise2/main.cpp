#include <iostream>
#include <cmath>
#include <algorithm>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;


//TEST 1: Verifica che tutti i marker siano associati correttamente

bool TestMarkers(const PolygonalMesh& mesh)
{
    // Stampa tutti i marker
    cout << "Marker registrati (punti):" << endl;

    for (const auto& [marker, list_id] : mesh.Cell0DMarkers)
    {
        cout << "Marker0D: " << marker << "  IDs = [";
        for (auto& id : list_id)
            cout << ' ' << id;
        cout << " ]" << endl;
    }

    cout << endl;
    cout << "Marker registrati (segmenti):" << endl;

    for (const auto& [marker, list_id] : mesh.Cell1DMarkers)
    {
        cout << "Marker1D: " << marker << "  IDs = [";
        for (const auto& id : list_id)
            cout << ' ' << id;
        cout << " ]" << endl;
    }
    
    cout << endl;    
    cout<<"TEST 1 SUPERATO: tutti i marker sono associati correttamente"<<endl;
    return true;
}

// Test 2: Verifica che ogni lato abbia lunghezza > 0

bool TestEdgesNonZero(const PolygonalMesh& mesh)
{
    for (unsigned int i = 0; i < mesh.NumCell2Ds; i++)  //itero sui poligoni della mesh
        for(unsigned int j = 0; j<mesh.Cell2DsEdges[i].size();j++) //itero sui lati del poligono
        {
            const vector<unsigned int>& lati = mesh.Cell2DsEdges[i];
            
            // ogni lato è segmento che congiunge due punti
            const int& punto1 = mesh.Cell1DsExtrema(0,lati[j]);  //origine 
            const int& punto2 = mesh.Cell1DsExtrema(1,lati[j]);  //fine lato

            const vector<double>& coordinate1 = {mesh.Cell0DsCoordinates(0,punto1),mesh.Cell0DsCoordinates(1,punto1)};
            const vector<double>& coordinate2 = {mesh.Cell0DsCoordinates(0,punto2),mesh.Cell0DsCoordinates(1,punto2)};
            
            double dx = coordinate2[0] - coordinate1[0];
            double dy = coordinate2[1] - coordinate1[1];
            double distanza = sqrt(dx*dx + dy*dy);

            if (distanza < 1e-16) 
            {
                cout<<"Errore nel poligono con ID "<< i <<", il lato "<< lati[j]<<" ha lunghezza 0"<<endl;
                return false;
            }
        }
    cout<<"TEST 2 SUPERATO: Verifica che ogni lato abbia lunghezza > 0"<<endl;
    return true;
}

// Test 3: Verifica che ogni poligono abbia area > 0

bool TestNonZeroArea(const PolygonalMesh& mesh)
{
    for (unsigned int i = 0; i < mesh.NumCell2Ds; i++) //itero sui poligoni della mesh
    {
        const auto& Idvertici = mesh.Cell2DsVertices[i];

        if (Idvertici.size() < 3)
            return false;

        double area = 0.0;

        for (size_t j = 0; j < Idvertici.size(); j++)
        {
            // id punti che formano l'i-esimo poligono
            unsigned int id1 = Idvertici[j];
            unsigned int id2 = Idvertici[(j + 1) % Idvertici.size()]; //j+1 modulo n

            const double& x1 = mesh.Cell0DsCoordinates(0,id1);
			const double& y1 = mesh.Cell0DsCoordinates(1,id1);
			const double& x2 = mesh.Cell0DsCoordinates(0,id2);
			const double& y2 = mesh.Cell0DsCoordinates(1,id2);
			
			area += x1*y2-x2*y1;
		}
		
		area = 0.5*abs(area);
	
		if(area < 1e-12) //tolleranza
		{
			cout<<"Errore nel poligono con ID "<<i<<": area nulla"<<endl;
			return false;
		}
	}
    cout<<"TEST 3 SUPERATO: Verifica che ogni poligono abbia area > 0"<<endl;
    return true;
}


int main()
{
    PolygonalMesh mesh;


    // Importazione mesh da file CSV
    if (!ImportMesh(mesh)) {
        cerr << "file not found" << endl;
        return 1;
    }

    /// Per visualizzare online le mesh:
    /// 1. Convertire i file .inp in file .vtu con https://meshconverter.it/it
    /// 2. Caricare il file .vtu su https://kitware.github.io/glance/app/

    // Test 1 – Verifica marker
    if (!TestMarkers(mesh))
    {
        cerr << "Errore: marker non validi." << endl;
        return 2;
    }

    // Test 2 – Verifica che ogni lato abbia lunghezza > 0
    if (!TestEdgesNonZero(mesh))
    {
        cerr << "Errore: esistono spigoli di lunghezza nulla." << endl;
        return 3;
    }

    // Test 3 – Verifica che ogni poligono abbia area > 0
    if (!TestNonZeroArea(mesh))
    {
        cerr << "Errore: esistono poligoni con area nulla." << endl;
        return 4;
    }

    // Test 4: Esportazione per verifica manuale
    Gedim::UCDUtilities utilities;
    {
        vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

        cell0Ds_properties[0].Label = "Marker";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);
        for(const auto &m : mesh.Cell0DMarkers)
            for(const unsigned int id: m.second)
                cell0Ds_marker.at(id) = m.first;

        cell0Ds_properties[0].Data = cell0Ds_marker.data();

        utilities.ExportPoints("./Cell0Ds.inp",
                               mesh.Cell0DsCoordinates,
                               cell0Ds_properties);
    }

    {

        vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].Label = "Marker";
        cell1Ds_properties[0].UnitLabel = "-";
        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);
        for(const auto &m : mesh.Cell1DMarkers)
            for(const unsigned int id: m.second)
                cell1Ds_marker.at(id) = m.first;

        cell1Ds_properties[0].Data = cell1Ds_marker.data();

        utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {},
                                 cell1Ds_properties);
    }

    cout << "Mesh caricata e validata con successo!" << endl;
    return 0;
}

