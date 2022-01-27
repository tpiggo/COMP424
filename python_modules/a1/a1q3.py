import math
import random
import itertools


class City:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return other.x == self.x and other.y == self.y


def city_exists(cities, x, y):
    for city in cities:
        if city.x == x and city.y == y:
            return True
    return False


def recursive_permutations(element_to_perm, gen_perm, current_perm):
    if element_to_perm:
        for element in element_to_perm:
            next_perm = current_perm + [element]
            remaining_perms = [item for item in element_to_perm if item != element]
            recursive_permutations(remaining_perms, gen_perm, next_perm)
    else:
        gen_perm.append(current_perm)
    return gen_perm


class TSPInstance:
    def __init__(self, cities):
        self._cities = cities
        self.city_permutations = []
        self.tours = []
        self.mean = 0
        self.max = 0
        self.path_sum = 0
        self.min = 0

    def brute_force_solve(self):
        self.city_permutations = itertools.permutations(self._cities, len(self._cities))
        for cities_order in self.city_permutations:
            number_cities = len(cities_order)
            tour_distance = 0
            for index, city in enumerate(cities_order):
                if number_cities == index + 1:
                    next_city = cities_order[0]
                else:
                    next_city = cities_order[index + 1]
                x_distance = next_city.x - city.x
                y_distance = next_city.y - city.y
                tour_distance += math.sqrt(math.pow(x_distance, 2) + math.pow(y_distance, 2))
            self.tours.append(tour_distance)
        self.calculate_reporting()

    def calculate_reporting(self):
        self.max = self.tours[0]
        self.min = self.tours[0]
        self.path_sum = self.tours[0]
        i = 1
        number_tours = len(self.tours)
        while i < number_tours:
            if self.min > self.tours[i]:
                self.min = self.tours[i]
            elif self.max < self.tours[i]:
                self.max = self.tours[i]
            self.path_sum += self.tours[i]
            i += 1
        self.mean = self.path_sum/number_tours

    @staticmethod
    def generate_tsp_instance(n=7):
        """
        Factory method for generating TSP instance
        :param n:
        :return:
        """
        cities = []
        _x = _y = 0
        for _ in range(0, n):
            is_gen = False
            while not is_gen:
                _x = random.random()
                _y = random.random()
                # is slow :(
                is_gen = not city_exists(cities, _x, _y)
            cities.append(City(_x, _y))
        return TSPInstance(cities)


class AggregateTSPInstance:
    def __init__(self, tsp_instances=None):
        self._tsp_instances = tsp_instances or []
        self.mean = 0
        self.max = 0
        self.min = 0
        self.path_sum = 0

    def brute_force_all(self):
        for tsp_inst in self._tsp_instances:
            tsp_inst.brute_force_solve()
        self.max = self._tsp_instances[0].min
        self.min = self._tsp_instances[0].min
        path_sum = self._tsp_instances[0].min
        num_instance = len(self._tsp_instances)
        for i in range(1, num_instance):
            if self.min > self._tsp_instances[i].min:
                self.min = self._tsp_instances[i].min
            elif self.max < self._tsp_instances[i].min:
                self.max = self._tsp_instances[i].min
            path_sum += self._tsp_instances[i].min
        self.mean = path_sum/num_instance
        # get std of optimal tour length
        self.std = 0


    def add_tsp_instance(self, tsp: TSPInstance):
        self._tsp_instances.append(tsp)

    @staticmethod
    def create_aggregate_tsp_instance(n=100):
        agg_tsp = AggregateTSPInstance()
        for _ in range(0, n):
            agg_tsp.add_tsp_instance(TSPInstance.generate_tsp_instance())
        return agg_tsp


if __name__ == '__main__':
    list_items = [1,2,3,4,5,6,7]
    # out = recursive_brute_force([1,2,3])
    out = recursive_permutations(list_items, [], [])
    print(len(out))
    #
    # agg_tsp_instance = AggregateTSPInstance.create_aggregate_tsp_instance()
    # agg_tsp_instance.brute_force_all()
    # tsp = TSPInstance.generate_tsp_instance()
    # tsp.brute_force_solve()