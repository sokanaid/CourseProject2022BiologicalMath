import pandas as pd
import os


def read_csv_file(directory_path, name):
    file_path = os.path.join(directory_path, name + ".csv")
    # print("read file ", file_path)
    return pd.read_csv(file_path)


def read_description(directory_path, file_name="description.txt", parse_params=False):
    file_path = os.path.join(directory_path, file_name)
    # print("read description from file  ", file_path)
    with open(file_path, 'r') as file:
        if not parse_params:
            return file.read()
        result = {}
        for line in file:
            name, value = line.split('=')
            value = value.strip(' \t\n\r')
            if value == 'TRUE':
                value = True
            elif value == 'FALSE':
                value = False
            else:
                try:
                    value = float(value)
                except:
                    try:
                        value = int(value)
                    except:
                        ...
            result[name] = value
        return result


def read_files(directory_path):
    data = {}
    # строковое описание параметров
    data["description_str"] = read_description(directory_path)
    # описание параметров в виде словаря
    data["params"] = read_description(directory_path, parse_params=True)
    # численность популяции
    data["population"] = read_csv_file(directory_path, "population")
    data["population"].drop(data["population"].tail(1).index, inplace=True)
    data["population"].rename(columns={"Unnamed: 0": 'epochs'}, inplace=True)
    # результаты последовательных сглаживаний численности популяции
    name = "exp_pop10"
    data[name] = read_csv_file(directory_path, name)
    data[name].rename(columns={"Unnamed: 0": 'epochs'}, inplace=True)
    return data


# читаем все результаты симуляций из папок
def read_all_simulations(directories_paths):
    data = []
    for directory_path in directories_paths:
        for simulation_diractory_name in os.listdir(directory_path):
            simulation_diractory_path = os.path.join(directory_path, simulation_diractory_name)
            if os.path.isdir(simulation_diractory_path):
                data.append(read_files(simulation_diractory_path))
    return data


def add_sim(data_frame, sim):
    params = sim["params"]
    mean_pop = sim["population"].tail(100)['pop'].mean()
    d = params['d'] if 'd' in params else 0
    new_row = pd.Series(
        {'b': params['b'], 'd': d, 'death_r': params['death_r'], 'dd': params['dd'], 'sd_b': params['sd_b'],
         'sd_d': params['sd_d'],
         'area_length_x': params['area_length_x'], 'initial_pop': params['initial_pop'], 'plateau_pop': mean_pop})
    data_frame = pd.concat([data_frame, new_row.to_frame().T], ignore_index=True)
    return data_frame


def main():
    directories = ["../../simulations_tables/data_set"]
    data = read_all_simulations(directories)
    sim_frame = pd.DataFrame(
        columns=['b', 'd', 'death_r', 'dd', 'sd_b', 'sd_d', 'area_length_x', 'initial_pop', 'plateau_pop'])
    for sim in data:
        if not sim["params"]["auto_stop_at_plateau"]:
            continue
        sim_frame = add_sim(sim_frame, sim)
    file_name = "data_set1.csv"
    sim_frame.to_csv(file_name, sep='\t')


if __name__ == '__main__':
    main()
