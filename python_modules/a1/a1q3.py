import datetime
import math
import random
import itertools
import statistics as sts
from typing import List
import datetime as dt


def gen_random(max_num):
    i = int(random.random() * max_num)
    i2 = int(random.random() * max_num)
    while i == i2:
        i2 = int(random.random() * max_num)
    return i, i2


def check_in(list_of_items, element1, element2):
    for x, y in list_of_items:
        if (x == element1 and y == element2) or (x == element2 and y == element1):
            return True
    return False


def recursive_permutations(element_to_perm, gen_perm, current_perm):
    if len(element_to_perm) > 0:
        for element in element_to_perm:
            next_perm = current_perm + [element]
            remaining_perms = [item for item in element_to_perm if item != element]
            recursive_permutations(remaining_perms, gen_perm, next_perm)
    else:
        gen_perm.append(current_perm)
    return gen_perm


def find_random_path(elements, order):
    if len(elements) == 0:
        return order
    elif len(elements) == 1:
        order.append(elements[0])
        return order
    next_choice = random.random() * len(elements)
    next_choice = elements[int(next_choice)]
    elements = [element for element in elements if element != next_choice]
    order.append(next_choice)
    return find_random_path(elements, order)


def find_min_max_mean(input_list):
    lmin = min(input_list)
    lmax = max(input_list)
    mean = sts.mean(input_list)
    std = sts.stdev(input_list, mean)
    return lmin, lmax, mean, std


class City:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return other.x == self.x and other.y == self.y

    def __ne__(self, other):
        return other.x != self.x and other.y != self.y

    def distance_to_next(self, next_city):
        return math.sqrt(math.pow(next_city.y - self.y, 2) + math.pow(next_city.y - self.y, 2))


class BaseInstance:
    def __init__(self):
        self.min = 0
        self.max = 0
        self.mean = 0
        self.std = 0

    def print_stats(self):
        print('For the following run we got:')
        print(f'Min: {self.min}')
        print(f'Max: {self.max}')
        print(f'Mean: {self.mean}')
        print(f'Standard deviation: {self.std}')


class TSPInstance(BaseInstance):
    def __init__(self, cities):
        super(TSPInstance, self).__init__()
        self._cities = cities
        self.city_permutations = []
        self.tours = []
        self.random_path_length = 0
        self.random_path = None
        self.hill_climbed_path_length = 0
        self._start_city = self._cities[0]

    def get_num_cities(self):
        return len(self._cities)

    def _calculate_tour(self, cities_order):
        if self._start_city == cities_order[0]:
            slice_index = 0
        elif self._start_city == cities_order[-1]:
            slice_index = len(cities_order) - 1
        else:
            slice_index = 0
            for index, city in enumerate(cities_order):
                if city == self._start_city:
                    slice_index = index
                    break
        num_cities = len(cities_order)
        max_num = slice_index + num_cities
        tour_distance = 0
        for index in range(slice_index, max_num):
            city = cities_order[index % num_cities]
            next_city = cities_order[(index + 1) % num_cities]
            tour_distance = tour_distance + city.distance_to_next(next_city)
        return tour_distance

    def brute_force_solve(self):
        self.city_permutations = itertools.permutations(self._cities, len(self._cities))
        for cities_order in self.city_permutations:
            self.tours.append(self._calculate_tour(cities_order))
        self.calculate_reporting()

    def calculate_reporting(self):
        self.max = max(self.tours)
        self.min = min(self.tours)

    def find_random_path(self):
        self.random_path = find_random_path(self._cities, [])

    def find_random_path_length(self):
        # No reporting here, only returning one
        self.find_random_path()
        self.random_path_length = self._calculate_tour(self.random_path)

    def solve_brute_with_random(self):
        self.brute_force_solve()
        self.find_random_path_length()

    def _swap_city(self, swap_index1, swap_index2):
        path = []
        for index, city in enumerate(self.random_path):
            if index == swap_index1:
                path.append(self.random_path[swap_index2])
            elif index == swap_index2:
                path.append(self.random_path[swap_index1])
            else:
                path.append(city)
        return path, self._calculate_tour(path)

    def _compute_all_neighbours(self):
        neighbours = []
        num_cities = len(self.random_path)
        for index in range(0, num_cities - 1):
            for index1 in range(index + 1, num_cities):
                neighbours.append(self._swap_city(index, index1))
        return neighbours

    def _compute_neighbour(self, computed):
        max_combs = math.comb(len(self.random_path), 2)
        if len(computed) == max_combs:
            return None
        max_length = len(self._cities)
        i, i2 = gen_random(max_length)
        while check_in(computed, i, i2):
            i, i2 = gen_random(max_length)
        computed.append((i, i2))
        return self._swap_city(i, i2), computed

    def solve_hill_climbing(self):
        self.find_random_path_length()
        current_tour_length = self.random_path_length
        computed = []
        while True:
            result = self._compute_neighbour(computed)
            if result is None:
                # Gone through all possible paths and found nothing better!
                break
            path_info, computed = result
            if path_info[1] < current_tour_length:
                self.random_path = path_info[0]
                current_tour_length = path_info[1]
        self.hill_climbed_path_length = current_tour_length

    @staticmethod
    def generate_tsp_instance(n: int):
        """
        Factory method for generating TSP instance
        :param n:
        :return:
        """
        cities = []
        new_city = None
        for _ in range(0, n):
            is_gen = False
            while not is_gen:
                new_city = City(random.random(), random.random())
                is_gen = new_city not in cities
            cities.append(new_city)
        return TSPInstance(cities)


class AggregateTSPInstance(BaseInstance):
    def __init__(self, tsp_instances: List[TSPInstance]):
        super(AggregateTSPInstance, self).__init__()
        self._tsp_instances = tsp_instances
        self.r_std = 0
        self.r_min = 0
        self.r_max = 0
        self.r_mean = 0
        self.optimal_rand = 0
        self.optimal_hc = 0
        self.hc_mean = 0
        self.hc_min = 0
        self.hc_max = 0
        self.hc_std = 0

    def _find_data_points(self):
        result_tuple = find_min_max_mean([tsp.min for tsp in self._tsp_instances])
        self.min, self.max, self.mean, self.std = result_tuple

    def _find_data_points_rand(self):
        result_tuple = find_min_max_mean([tsp.random_path_length for tsp in self._tsp_instances])
        # get std of optimal tour length
        self.r_min, self.r_max, self.r_mean, self.r_std = result_tuple

    def _find_data_points_hc(self):
        result_tuple = find_min_max_mean([tsp.hill_climbed_path_length for tsp in self._tsp_instances])
        # get std of optimal tour length
        self.hc_min, self.hc_max, self.hc_mean, self.hc_std = result_tuple

    def find_optimal_vs_rand(self):
        self.optimal_rand = 0
        for inst in self._tsp_instances:
            if inst.min == inst.random_path_length:
                self.optimal_rand += 1

    def find_optimal_vs_hill_climbing(self):
        self.optimal_hc = 0
        for inst in self._tsp_instances:
            if inst.min == inst.hill_climbed_path_length:
                self.optimal_hc += 1

    def brute_force_all(self):
        if len(self._tsp_instances) == 0:
            raise ValueError("Missing instances to run this over!")
        for tsp_inst in self._tsp_instances:
            tsp_inst.brute_force_solve()
        self._find_data_points()

    def print_stats_rand(self, check_optimal=True):
        print('Random findings:')
        print(f'Min: {self.r_min}')
        print(f'Max: {self.r_max}')
        print(f'Mean: {self.r_mean}')
        print(f'Standard deviation: {self.r_std}')
        if check_optimal:
            print(f'Number of optimal: {self.optimal_rand}')

    def print_stats_hill_climbing(self, check_optimal=True):
        print('Hill Climbing findings:')
        print(f'Min: {self.hc_min}')
        print(f'Max: {self.hc_max}')
        print(f'Mean: {self.hc_mean}')
        print(f'Standard deviation: {self.hc_std}')
        if check_optimal:
            print(f'Number of optimal: {self.optimal_hc}')

    def random_paths(self):
        for tsp_inst in self._tsp_instances:
            tsp_inst.solve_brute_with_random()
        self._find_data_points()
        self._find_data_points_rand()
        self.find_optimal_vs_rand()
        
    def find_random_paths(self):
        for tsp_inst in self._tsp_instances:
            tsp_inst.find_random_path_length()
        self._find_data_points_rand()

    def solve_hill_climbing(self, print_done=False):
        print(f"Starting hill climbing for {self._tsp_instances[0].get_num_cities()} instances!")
        start = datetime.datetime.now()
        for index, tsp_inst in enumerate(self._tsp_instances):
            tsp_inst.solve_hill_climbing()
            if print_done:
                print('Done hill climbing for iteration {}'.format(index))
        self._find_data_points_hc()
        print("Done hill climbing")
        print(f"Took {datetime.datetime.now() - start}")

    def brute_and_hillclimbing(self):
        for tsp_inst in self._tsp_instances:
            tsp_inst.brute_force_solve()
            tsp_inst.solve_hill_climbing()
        self._find_data_points()
        self._find_data_points_hc()
        self.find_optimal_vs_hill_climbing()

    @staticmethod
    def create_aggregate_tsp_instance(n=100, num_cities=7):
        agg_tsp = []
        for _ in range(0, n):
            agg_tsp.append(TSPInstance.generate_tsp_instance(num_cities))
        return AggregateTSPInstance(agg_tsp)


def run_brute_force(n=7):
    agg_tsp_instance = AggregateTSPInstance.create_aggregate_tsp_instance(num_cities=n)
    agg_tsp_instance.brute_force_all()
    agg_tsp_instance.print_stats()


def run_random_vs_brute(n=7):
    agg_tsp_instance = AggregateTSPInstance.create_aggregate_tsp_instance(num_cities=n)
    agg_tsp_instance.random_paths()
    agg_tsp_instance.print_stats()
    agg_tsp_instance.print_stats_rand()


def run_random(n=100):
    agg_tsp_instance = AggregateTSPInstance.create_aggregate_tsp_instance(num_cities=n)
    agg_tsp_instance.find_random_paths()
    agg_tsp_instance.print_stats_rand(False)


def run_hill_climbing_vs_brute(n=7):
    agg_tsp_instance = AggregateTSPInstance.create_aggregate_tsp_instance(num_cities=n)
    agg_tsp_instance.brute_and_hillclimbing()
    agg_tsp_instance.print_stats()
    agg_tsp_instance.print_stats_hill_climbing()


def run_hill_climbing(n=100):
    agg_tsp_instance = AggregateTSPInstance.create_aggregate_tsp_instance(num_cities=n)
    agg_tsp_instance.solve_hill_climbing(True)
    agg_tsp_instance.print_stats_hill_climbing(False)


def main():
    # run_brute_force()
    # run_random_vs_brute()
    # run_hill_climbing_vs_brute()
    run_random()
    run_hill_climbing()


if __name__ == '__main__':
    main()