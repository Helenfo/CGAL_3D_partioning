//#define CGAL_USE_BASIC_VIEWER

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>
#include <list>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/convexity_check_3.h>
#include <fstream>
#include <iostream>
using std::vector;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;
typedef Polyhedron_3::HalfedgeDS             HalfedgeDS;
typedef typename HalfedgeDS::Vertex   Vertex;
typedef typename Vertex::Point Point;
typedef Polyhedron_3::Vertex_iterator        Vertex_iterator;
typedef Polyhedron_3::Facet_iterator        Facet_iterator;
typedef Polyhedron_3::Halfedge_around_facet_circulator HF_circulator;
typedef std::vector<int> Facet;
typedef std::vector<Point> VertexList;
typedef std::vector<Facet> FacetList;


struct MeshPolyhedron {
	VertexList vertices;
	FacetList facets;

	MeshPolyhedron() {}
	MeshPolyhedron(VertexList vs, FacetList fs) {
		this->vertices = vs;
		this->facets = fs;
	}
};

template <class HDS>
class BuildPolyhedron : public CGAL::Modifier_base<HDS> {
private:
	VertexList vertices;
	FacetList facets;
public:
	BuildPolyhedron(VertexList vertices, FacetList facets)
	{

		this->vertices.insert(std::end(this->vertices),
			std::begin(vertices), std::end(vertices));
		this->facets.insert(std::end(this->facets),
			std::begin(facets), std::end(facets));
	}

	void operator()(HDS& hds) {
		// Postcondition: hds is a valid polyhedral surface.
		CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
		int num_vertices = vertices.size();
		int num_facets = facets.size();
		B.begin_surface(num_vertices, num_facets);
		typedef typename HDS::Vertex   Vertex;
		typedef typename Vertex::Point Point;
		for (auto v : vertices)
			B.add_vertex(v);

		int k = 0;
		for (auto facet : facets) {
			k++;
			bool isfacet = B.test_facet(facet.begin(), facet.end());
			if (!isfacet) {
				std::cout << "obs, facet not well defined " << k << "\n.";
			}
			B.begin_facet();
			for (auto point : facet) {
				B.add_vertex_to_facet(point);
			}
			B.end_facet();
		}
		B.end_surface();

	}
};

Polyhedron_3 createPolyhedron(double vs[][3], int num_vertices, int fs[][8], int num_facets, std::vector<vector<int>> faces);
std::list<Polyhedron_3> decompose(Nef_polyhedron_3 N);
std::list<Polyhedron_3> decompose(Polyhedron_3 P);
void poly2mesh(Polyhedron_3 P, MeshPolyhedron& MP);
int find(VertexList vs, Point p);

int main() {
	// obstacle info
	int num_vertices = 16;
	int num_facets = 10;

	//double v[12][3] = { {0, 0, 0}, { 50,0,0 }, { 50,50,0 }, { 0,50,0 }, { 50,0,25 }, { 50,50,25 }, { 25,0,25 }, { 25,50,25 }, { 25,0,50 }, { 25,50,50 }, { 0,0,50 }, { 0,50,50 } };
	//int f[8][6] = { {1,11,9,7,5,2}, {6,8,10,12,4,3},{11,1,4,12},{5,6,3,2},{12,10,9,11},{3,4,1,2},{9,10,8,7},{8,6,5,7} };
	double v[16][3] = { {0,0,0}, {10,0,0}, {10,10,0}, {0,10,0}, {20,0,0},{30,0,0,}, {30,10,0}, {20,10,0}, {10,0,10}, {20,0,10}, {20,10,10}, {10,10,10}, {0,0,20}, {30,0,20}, {30,10,20}, {0,10,20} };
	int f[10][8] = { {16,15,14,13}, {14,15,7,6}, {16,13,1,4}, {9,10,11,12}, {5,6,7,8}, {1,2,3,4}, {9,12,3,2}, {11,10,5,8}, {15,16,4,3,12,11,8,7}, {13,14,6,5,10,9,2,1} };

	//double v[16][3] = {{ 0, 0, 0 }, { 30,0,0 }, { 30,10,0 }, { 0,10,0 }, { 10,0,10 }, { 20,0,10 }, { 20,10,10 }, { 10,10,10 }, { 10,0,20 }, { 20,0,20 }, { 20, 10, 20 }, { 10,10,20 }, { 30,0,30 }, { 30,10,30 }, { 0,0,30 }, { 0,10,30 }};
	//int f[10][8] = { {1,2,3,4}, {16,14,13,15}, {9,10,11,12}, {8,7,6,5}, {10,6,7,11},{9,12,8,5}, {13,14,3,2}, {16,15,1,4}, { 15,13,9,10,6,2,1,5 }, {14,16,4,3,11,12,8,7} };

	std::vector<vector<int>> faces = { {16,15,14,13}, {14,15,7,6}, {16,13,1,4}, {9,10,11,12}, {5,6,7,8}, {1,2,3,4}, {9,12,3,2}, {11,10,5,8}, {15,16,4,3,12,11,8,7}, {13,14,6,5,10,9,2,1} };
	//std::vector<vector<int>> faces = { {1, 2, 3, 4}, { 16,14,13,15 }, { 9,10,11,12 }, { 8,7,6,5 }, { 10,6,7,11 }, { 9,12,8,5 }, { 13,14,3,2 }, { 16,15,1,4 }, { 15,13,9,10,6,2,1,5 }, { 14,16,4,3,11,12,8,7 } };

	Polyhedron_3 Pobs = createPolyhedron(v, num_vertices, f, num_facets, faces);

	Nef_polyhedron_3 Nobs(Pobs);

	std::list<Polyhedron_3> convex_parts = decompose(Nobs);
	//Polyhedron_3 P1 = convex_parts;

	std::list<Polyhedron_3>::iterator it;

	//Polyhedron_3 P2;
	std::ofstream Hypnos_FILE;
	std::string TEXTO = "Escrevendo em arquivo de texto";
	Hypnos_FILE.open("C:\\Users\\helen\\Desktop\\CGAL_ok\\CGAL-5.3-examples\\CGAL-5.3\\examples\\Convex_decomposition_3\\build\\Debug\\Output_partioning_final.txt", std::ios::app);
	if (Hypnos_FILE.is_open())
	{
		std::cout << "Arquivo de texto aberto com sucesso!" << std::endl;

		for (it = convex_parts.begin(); it != convex_parts.end(); it++)
		{
			Polyhedron_3 P2 = *it;
			MeshPolyhedron MP;
			poly2mesh(P2, MP);

			std::cout << "polyhedro all facets\n";
			Hypnos_FILE << "\n[";
			for (auto f : MP.facets) {
				std::cout << "nova face\n";
				Hypnos_FILE << ";";
				for (auto i : f) {
					std::cout << MP.vertices[i] << std::endl;
					Hypnos_FILE << MP.vertices[i] << ";";
				}
				std::cout << "\n";
				//Hypnos_FILE << "\n";
			}
			std::cout << "\n";
			Hypnos_FILE << "]\n";
		}

	}
	else
		std::cout << "Erro ao abrir arquivo de texto.";

	Hypnos_FILE.close();
	

}

void poly2mesh(Polyhedron_3 P, MeshPolyhedron& MP)
{
	std::size_t i = 0;
	for (Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v)
		MP.vertices.push_back(v->point());

	for (Facet_iterator f = P.facets_begin(); f != P.facets_end(); ++f) {
		i++;
		HF_circulator hinit = f->facet_begin();
		HF_circulator h = hinit;
		Facet facet;
		do {
			Point p = h->vertex()->point();
			int index = find(MP.vertices, p);
			facet.push_back(index);
			h++;
		} while (h != hinit);
		MP.facets.push_back(facet);
	}
}

int find(VertexList vs, Point p)
{
	int index = 0;
	for (auto v : vs) {
		if (p[0] == v[0] && p[1] == v[1] && p[2] == v[2])
			return index;
		index++;
	}
	return -1;
}

Polyhedron_3 createPolyhedron(double vs[][3], int num_vertices, int fs[][8], int num_facets, std::vector<vector<int>> faces)
{
	VertexList vertices;
	for (int i = 0; i < num_vertices; ++i)
		vertices.push_back(Point(vs[i][0], vs[i][1], vs[i][2]));
	FacetList facets;
	for (int i = 0; i < num_facets; ++i) {
		facets.push_back(Facet());
		for (int j = 0; j < faces[i].size(); ++j)
			facets[i].push_back(faces[i][j] - 1);
	}

	Polyhedron_3 P;
	BuildPolyhedron<HalfedgeDS> poly(vertices, facets);
	P.delegate(poly);

	return P;
}

std::list<Polyhedron_3> decompose(Polyhedron_3 P)
{
	Nef_polyhedron_3 N(P);
	return decompose(N);
}

std::list<Polyhedron_3> decompose(Nef_polyhedron_3 N)
{
	CGAL::convex_decomposition_3(N);
	std::list<Polyhedron_3> convex_parts;

	// the first volume is the outer volume, which is 
	// ignored in the decomposition
	Volume_const_iterator ci = ++N.volumes_begin();
	for (; ci != N.volumes_end(); ++ci) {
		if (ci->mark()) {
			Polyhedron_3 P;
			N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
			convex_parts.push_back(P);
		}
	}
	std::cout << "decomposition into " << convex_parts.size() << " convex parts " << std::endl;

	return convex_parts;
}