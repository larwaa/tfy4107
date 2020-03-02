import statistics

CUTOFF = 1


def csv_to_list(csv):
    result = []
    with open(csv, "r") as file:
        file.readline()
        tabs = file.readline().split(",")
        result.append(tabs)
        for line in file.readlines()[2:]:
            split_line = line.split(",")
            time = float(split_line[0].strip())
            x = float(split_line[1].strip())
            y = float(split_line[2].strip())
            try:
                velocity = float(split_line[3].strip())
                result.append([time, x, y, velocity])
            except ValueError:
                pass
    return result


def list_to_avg_velocity(values):
    velocities = []
    for value in values[1:]:
        velocity = value[3]
        try:
            curr_average = statistics.mean(velocities)
            if (curr_average * (1 - CUTOFF)) < velocity < (curr_average * (1 + CUTOFF)):
                velocities.append(velocity)
        except statistics.StatisticsError:
            velocities.append(velocity)
    return statistics.mean(velocities)


print("#" * 5 + "BALL" + "#"*5)
for i in range(1, 11):
    values = csv_to_list(f"b{i}-v")
    print(list_to_avg_velocity(values))

print()
print("#" * 5 + "CIRCLE" + "#"*5)
for i in range(1, 11):
    values = csv_to_list(f"c{i}-v")
    print(list_to_avg_velocity(values))



