from karraynav import *

if __name__ == "__main__":
    tests = kArrayTests()
    tests.tests_general()
    tests.tests_vector()
    tests.tests_matrix()

    tests = kArrayNavTests()
    tests.do_tests()
