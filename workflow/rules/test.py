print("Hello world")
import os

reference_panel_dirpath = "cesta/"
_SUFFIX = ".fa"
_REGEX = ".*"
location_format = os.path.join(reference_panel_dirpath, f"{{name, {_REGEX}}}{_SUFFIX}")
print(location_format)

location_format = f"{reference_panel_dirpath}/{{name, {_REGEX}}}{_SUFFIX}"
print(location_format)
