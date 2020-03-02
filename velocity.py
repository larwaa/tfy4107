import numpy as np
import matplotlib as plt

CUTOFF = 0.5


def csv_to_list(csv):
    result = []
    with open(csv, "r") as file:
        tabs = file.readline().split(",")
        result.append(tabs)
        for line in file[1:]:
            split_line = line.split(",")
            time = split_line[0]
            x = split_line[1]
            y = split_line[2]
            velocity = split_line[3]
            result.append(time, x, y, velocity)
    return result


def list_to_avg_velocity(values):
    velocities = []
    for value in values[1:]:
        velocity = value[3]
        curr_average = avg(velocities)
        if curr_average == 0 or ((curr_average * (1 - CUTOFF)) < velocities < (curr_average * (1 + CUTOFF))):
            velocities.append(velocity)
    return avg(velocities)