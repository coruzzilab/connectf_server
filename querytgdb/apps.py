import matplotlib
from django.apps import AppConfig

matplotlib.use('SVG')

import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ["DejaVu Sans"]


class QuerytgdbConfig(AppConfig):
    name = 'querytgdb'
