from brian2 import *

def plot_variables_with_axis(plotting_list, sv_mon):
    count = len(plotting_list)
    fig, ax = plt.subplots(count, 1,sharex=True)
    ax[-1].set_xlabel("X-axis Label")
    a = 0
    for element in plotting_list:
        print(element['variable'])
        ax[a].plot(getattr(sv_mon, element['variable'])[:].T)
        ax[a].set_ylabel(element['axis'])
        a+=1