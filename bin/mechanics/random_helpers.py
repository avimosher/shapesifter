import random
import string

random.seed(812321)

def random_name():
    return ''.join(random.choice(string.ascii_lowercase+string.digits) for i in range(12))
def random_vector(low,high):
    return [random.uniform(l,h) for l,h in zip(low,high)]

