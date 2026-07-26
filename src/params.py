import yaml

def load_params(path="config/parameters.yml"):
    with open(path) as f:
        params = yaml.safe_load(f)

    half_width_m = params["domain"]["width_km"] * 1000 / 2
    half_height_m = params["domain"]["height_km"] * 1000 / 2

    params["domain"]["xmin"] = -half_width_m
    params["domain"]["xmax"] = half_width_m
    params["domain"]["ymin"] = -half_height_m
    params["domain"]["ymax"] = half_height_m
    params["domain"]["corners"] = [
        (-half_width_m, -half_height_m),
        (half_width_m, -half_height_m),
        (half_width_m, half_height_m),
        (-half_width_m, half_height_m),
    ]
    return params