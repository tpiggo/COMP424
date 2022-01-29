import math
import random
import itertools
from typing import List


class City:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return other.x == self.x and other.y == self.y

    def __ne__(self, other):
        return other.x != self.x and other.y != self.y

    def distance_to_next(self, next):
        x_distance = next.x - self.x
        y_distance = next.y - self.y
        return math.sqrt(math.pow(x_distance, 2) + math.pow(y_distance, 2))


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

class BaseInstance:
    def __init__(self):
        self.min = 0
        self.max = 0
        self.mean = 0
        self.std = 0
        self.path_sum = 0

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

    @staticmethod
    def _calculate_tour(cities_order):
        number_cities = len(cities_order)
        tour_distance = 0
        for index, city in enumerate(cities_order):
            if number_cities == index + 1:
                next_city = cities_order[0]
            else:
                next_city = cities_order[index + 1]
            tour_distance += city.distance_to_next(next_city)
        return tour_distance

    def brute_force_solve(self):
        self.city_permutations = itertools.permutations(self._cities, len(self._cities))
        for cities_order in self.city_permutations:
            self.tours.append(TSPInstance._calculate_tour(cities_order))
        self.calculate_reporting()

    def calculate_reporting(self):
        self.max = self.tours[0]
        self.min = self.tours[0]
        for tour in self.tours:
            if self.min > tour:
                self.min = tour
            elif self.max < tour:
                self.max = tour

    def find_random_path(self):
        self.random_path = find_random_path(self._cities, [])

    def find_random_path_length(self):
        # No reporting here, only returning one
        self.find_random_path()
        self.random_path_length = TSPInstance._calculate_tour(self.random_path)

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
        return path, len(path)

    def solve_hill_climbing(self):
        self.find_random_path()
        current_tour_length = TSPInstance._calculate_tour(self.random_path)
        # randomly choose a city
        first_index = int(random.random() * len(self.random_path))
        while True:
            next_index = first_index + 1 if first_index < len(self.random_path) - 1 else 0
            prev_index = first_index - 1 if first_index > 0 else len(self.random_path) - 1
            if next_index == len(self.random_path):
                print("GOING TO GO OUT OF BOUNDS!")
            npath, nlength = self._swap_city(first_index, next_index)
            ppath, plength = self._swap_city(first_index, prev_index)
            if (plength == nlength or plength > nlength) and nlength < current_tour_length:
                self.random_path = npath
                current_tour_length = TSPInstance._calculate_tour(self.random_path)
                first_index = next_index
            elif plength < current_tour_length:
                self.random_path = ppath
                current_tour_length = TSPInstance._calculate_tour(self.random_path)
                first_index = prev_index
            else:
                # Neither neighbour is optimal!
                break
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
                # is slow :(
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

    def _calculate_standard_deviation(self):
        total = 0.0
        for inst in self._tsp_instances:
            total = total + math.pow(inst.min - self.mean, 2)
        self.std = math.sqrt(total/len(self._tsp_instances))

    def _find_data_points(self):
        path_sum = 0
        for index, instance in enumerate(self._tsp_instances):
            if index == 0:
                self.max = instance.min
                self.min = instance.min
                path_sum = instance.min
            elif self.min > instance.min:
                self.min = instance.min
            elif self.max < instance.min:
                self.max = instance.min
            path_sum += instance.min
        self.mean = path_sum / len(self._tsp_instances)
        # get std of optimal tour length

    def _find_data_points_rand(self):
        path_sum = 0
        for index, instance in enumerate(self._tsp_instances):
            if index == 0:
                self.r_max = instance.random_path_length
                self.r_min = instance.random_path_length
                path_sum = instance.random_path_length
            elif self.r_min > instance.random_path_length:
                self.r_min = instance.random_path_length
            elif self.r_max < instance.random_path_length:
                self.r_max = instance.random_path_length
            path_sum += instance.random_path_length
        self.r_mean = path_sum / len(self._tsp_instances)
        # get std of optimal tour length

    def _find_data_points_hc(self):
        path_sum = 0
        for index, instance in enumerate(self._tsp_instances):
            if index == 0:
                self.hc_max = instance.hill_climbed_path_length
                self.hc_min = instance.hill_climbed_path_length
                path_sum = instance.hill_climbed_path_length
            elif self.hc_min > instance.hill_climbed_path_length:
                self.hc_min = instance.hill_climbed_path_length
            elif self.hc_max < instance.hill_climbed_path_length:
                self.hc_max = instance.hill_climbed_path_length
            path_sum += instance.hill_climbed_path_length
        self.hc_mean = path_sum / len(self._tsp_instances)
        # get std of optimal tour length

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

    def _calculate_standard_deviation_hc(self):
        total = 0.0
        for inst in self._tsp_instances:
            total = total + math.pow(inst.hill_climbed_path_length - self.hc_mean, 2)
        self.hc_std = math.sqrt(total / len(self._tsp_instances))

    def _calculate_standard_deviation_random(self):
        total = 0.0
        for inst in self._tsp_instances:
            total = total + math.pow(inst.random_path_length - self.r_mean, 2)
        self.r_std = math.sqrt(total/len(self._tsp_instances))

    def brute_force_all(self):
        if len(self._tsp_instances) == 0:
            raise ValueError("Missing instances to run this over!")
        for tsp_inst in self._tsp_instances:
            tsp_inst.brute_force_solve()
        self._find_data_points()
        # get std of optimal tour length
        self._calculate_standard_deviation()

    def print_stats(self, print_random=False, print_hc=False):
        super(AggregateTSPInstance, self).print_stats()
        if print_random:
            print('Random findings:')
            print(f'Min: {self.r_min}')
            print(f'Max: {self.r_max}')
            print(f'Mean: {self.r_mean}')
            print(f'Standard deviation: {self.r_std}')
            print(f'Number of optimal: {self.optimal_rand}')
        if print_hc:
            print('Hill Climbing findings:')
            print(f'Min: {self.hc_min}')
            print(f'Max: {self.hc_max}')
            print(f'Mean: {self.hc_mean}')
            print(f'Standard deviation: {self.hc_std}')
            print(f'Number of optimal: {self.optimal_hc}')

    def random_paths(self):
        self.path_sum = 0
        for tsp_inst in self._tsp_instances:
            tsp_inst.solve_brute_with_random()
        self._find_data_points()
        self._calculate_standard_deviation()
        self._find_data_points_rand()
        self._calculate_standard_deviation_random()
        self.find_optimal_vs_rand()

    def brute_and_hillclimbing(self):
        self.path_sum = 0
        for tsp_inst in self._tsp_instances:
            tsp_inst.solve_hill_climbing()
            tsp_inst.brute_force_solve()
        self._find_data_points()
        self._calculate_standard_deviation()
        self._find_data_points_hc()
        self._calculate_standard_deviation_hc()
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


def run_random_path(n=7):
    agg_tsp_instance = AggregateTSPInstance.create_aggregate_tsp_instance(num_cities=n)
    agg_tsp_instance.random_paths()
    agg_tsp_instance.print_stats(True)


def run_hill_climbing_single(n=7):
    tsp = TSPInstance.generate_tsp_instance(n)
    tsp.solve_hill_climbing()
    tsp.brute_force_solve()
    tsp.print_stats()
    print(f"Length of the hill climbed path {tsp.hill_climbed_path_length}")

def run_hill_climbing(n=7):
    agg_tsp_instance = AggregateTSPInstance.create_aggregate_tsp_instance(num_cities=n)
    agg_tsp_instance.brute_and_hillclimbing()
    agg_tsp_instance.print_stats(print_hc=True)

def main():
    # run_brute_force()
    # run_random_path()
    # run_hill_climbing_single()
    run_hill_climbing()


def test():
    tsp_inst = TSPInstance.generate_tsp_instance(4)
    distance = tsp_inst.find_random_path()
    print(f'random path distance is {distance}')


if __name__ == '__main__':
    main()
    # test()