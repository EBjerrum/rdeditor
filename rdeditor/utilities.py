import logging


def validate_rgb(rgb: tuple) -> None:
    if not isinstance(rgb, tuple) or len(rgb) != 3:
        logging.error(f"Input must be a tuple of length 3 representing an RGB color, was {rgb}")
        return False
        # raise ValueError("Input must be a tuple of length 3 representing an RGB color")

    for value in rgb:
        if not (isinstance(value, float) or isinstance(value, int)) or value < 0 or value > 1:
            logging.error(f"RGB values must be floats or integers between 0 and 1, was {rgb}")
            return False
            # raise ValueError("RGB values must be floats or integers between 0 and 1")

    return True
