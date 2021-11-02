#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

using namespace std;

class Edge
{
public:
    float a;
    float b;
    float weight;
    float pheromone;

    Edge(float a = 0, float b = 0, float weight = 0, float initial_pheromone = 0) {
        this->a = a;
        this->b = b;
        this->weight = weight;
        this->pheromone = initial_pheromone;
    }

    Edge& operator=(const Edge& other) {
        a = other.a;
        b = other.b;
        weight = other.weight;
        pheromone = other.pheromone;
        return *this;
    }
};



class Ant
{
public:
    float alpha;
    float beta;
    int num_nodes;
    vector<vector<Edge>> edges;
    vector<int> tour;
    float distance;

    Ant(float alpha, float beta, int num_nodes, vector<vector<Edge>> edges) {
        this->alpha = alpha;
        this->beta = beta;
        this->num_nodes = num_nodes;
        this->edges = edges;
    }

    int _select_node() {
        float roulette_wheel = 0.0;
        float heuristic_total = 0.0;
        
        vector<int> unvisited_nodes;
        for(int i = 0; i < num_nodes; i++) {
            if(find(tour.begin(), tour.end(), i) == tour.end())
                unvisited_nodes.push_back(i);
        }

        for(auto& unvisited_node : unvisited_nodes) {
            try {
                heuristic_total += edges[*(tour.end()-1)][unvisited_node].weight;
            } catch(const std::exception& e) {
                cout << e.what() << endl;
            }
        }

        for(auto& unvisited_node : unvisited_nodes) {
            roulette_wheel += pow(edges[*(tour.end()-1)][unvisited_node].pheromone, alpha) *
                pow((heuristic_total / edges[*(tour.end()-1)][unvisited_node].weight), beta);
        }

        auto random_value = rand() % (int)roulette_wheel;

        float wheel_position = 0;
        for(auto& unvisited_node : unvisited_nodes) {
            wheel_position += pow(edges[*(tour.end()-1)][unvisited_node].pheromone, alpha) *
                pow((heuristic_total / edges[*(tour.end()-1)][unvisited_node].weight), beta);
            if(wheel_position >= random_value) {
                return unvisited_node;
            }
        }
    }

    vector<int> find_tour() {
        tour.clear();
        tour.push_back(rand() % num_nodes);
        while(tour.size() < num_nodes) {
            auto node = _select_node();
            tour.push_back(node);
        }
        return tour;
    }

    float get_distance() {
        distance = 0;
        for (int i = 0; i < num_nodes; i++) {
            distance += edges[tour[i]][tour[(i + 1) % num_nodes]].weight;
        }
        return distance;
    }
};

class SolveTSPUsingACO
{
public:
    int colony_size;
    float pheromone_deposit_weight;
    int steps;
    float rho;
    vector<pair<float, float>> nodes;
    vector<vector<Edge>> edges;
    int num_nodes;
    vector<Ant> ants;
    vector<int> global_best_tour;
    float global_best_distance;

    vector<float> distances_data;

    SolveTSPUsingACO(vector<pair<float, float>> nodes, int colony_size = 10, float alpha = 1.0, float beta = 3.0, float rho = 0.1, float pheromone_deposit_weight = 1.0, float initial_pheromone=1.0, int steps = 100) {
        this->nodes = nodes;
        this->colony_size = colony_size;
        this->steps = steps;
        this->rho = rho;
        this->pheromone_deposit_weight = pheromone_deposit_weight;
        this->num_nodes = nodes.size();

        edges = vector<vector<Edge>>(num_nodes, vector<Edge>(num_nodes, Edge()));
        for (int i = 0; i < num_nodes; i++){
            for(int j = i +1; j< num_nodes; j++) {
                edges[i][j] = edges[j][i] = Edge(i, j, sqrt(
                    pow(nodes[i].first - nodes[j].first, 2.0) + pow(nodes[i].second - nodes[j].second, 2.0)),
                                                                initial_pheromone);
            }
        }

        ants = vector<Ant>(colony_size, Ant(alpha, beta, num_nodes, edges));
        global_best_distance = numeric_limits<float>::max();
    }

    void _add_pheromone(vector<int> tour, float distance, float weight=1.0) {
        float pheromone_to_add = pheromone_deposit_weight / distance;
        for(int i = 0; i < num_nodes; i++) {
            try{
                edges[tour[i]][tour[(i + 1) % num_nodes]].pheromone += weight * pheromone_to_add;
            } catch(const std::exception& e) {
                cout << e.what() << endl;
                cout << tour[i] << endl;
                cout << tour[(i + 1) % num_nodes] << endl;
                cout << edges.size() << endl;
            }
            
        }
    }

    void _acs(float lmin) {
        for(int i = 0;global_best_distance > lmin && i < steps; i++) {
            for(auto& ant : ants) {
                auto tour = ant.find_tour();
                auto distance = ant.get_distance();
                _add_pheromone(tour, distance);

                if(ant.distance < global_best_distance) {
                    global_best_tour = ant.tour;
                    global_best_distance = ant.distance;
                }
            }
            for(int j = 0; j < num_nodes; j++) {
                for(int k = j + 1; k < num_nodes; k++) {
                    edges[j][k].pheromone *= (1.0 - rho);
                }
            }

            distances_data.push_back(global_best_distance);
        }
    }

    void run(float lmin) {
        auto begin = chrono::steady_clock::now();
        cout << "Run ACS over TSP" << endl;
        _acs(lmin);
        cout << "End" << endl;
        cout << "Time: ";
        cout << (float)chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - begin).count() / 1000.0 << " ms" << endl;
        for(auto item : global_best_tour) {
            cout << item << " ";
        }
        cout << endl;
    }

    void plot() {
        float line_width = 1;
        float point_radius = sqrt(2);
        int annotation_size = 8;
        int dpi = 120;
        bool save = true;

        vector<float> xs;
        vector<float> ys;
        for(int i : global_best_tour) {
            xs.push_back(nodes[i].first);
            ys.push_back(nodes[i].second);
        }
        xs.push_back(xs[0]);
        ys.push_back(ys[0]);

        plt::plot(xs, ys);
        plt::scatter(xs, ys, 3.14 * pow(point_radius,2));
        plt::title("ACS over TSP");
        plt::show();

        plt::plot(distances_data);
        plt::title("Min distance during time");
        plt::show();
    }
};

class SolveTSPUsingHungryAlgo
{
    vector<pair<float, float>> nodes;
    vector<pair<float, float>> tour;
public:
    SolveTSPUsingHungryAlgo(vector<pair<float, float>> nodes) {
        this->nodes = nodes;
    }

    // the function returns 
    float run() {
        vector<pair<float, float>> visited;

        auto first = nodes[rand() % nodes.size()];
        tour.push_back(first);
        nodes.erase(find(nodes.begin(), nodes.end(), first));

        float Lmin = 0;
        while(nodes.size()) {
            tour.push_back(get_next_min());
            Lmin += length_between(tour[tour.size() - 1], tour[tour.size() - 2]);
        }

        Lmin += length_between(tour[0], tour[tour.size() - 1]);
        
        return Lmin;
    }

private:
    pair<float, float> get_next_min(){
        auto last = tour[tour.size()-1];
        
        pair<float, float> min = nodes[0];
        for(auto n : nodes) {
            auto prev = length_between(last, min);
            auto next = length_between(last, n);
            if(next < prev) {
                min = n;
            }
        }

        nodes.erase(find(nodes.begin(), nodes.end(), min));

        return min;
    }

    float length_between(pair<float, float> n1, pair<float, float> n2) {
        return sqrt(pow(n1.first - n2.first, 2) + pow(n1.second - n2.second, 2));
    }
};

int main() {
    srand(time(0));

    int colony_size = 35;
    int steps = 100;
    int num_nodes = 150;
    vector<pair<float, float>> nodes;
    for(int i = 0; i < num_nodes; i++) {
        int x = (rand() % 800) - 400;
        int y = (rand() % 800) - 400;
        nodes.push_back(make_pair(float(x), float(y)));
    }

    SolveTSPUsingHungryAlgo ha(nodes);
    auto lmin = ha.run();
    cout << "lmin: " << lmin << endl;

    SolveTSPUsingACO aco(nodes, colony_size, 1, 3, 0.4, 1.0, 1.0, steps);
    aco.run(lmin);
    cout << "shortest: " << aco.global_best_distance << endl;
    aco.plot();
}
