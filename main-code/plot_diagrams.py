from matplotlib import pyplot as plt


def general_diagram(inserts_amount, replaces_amount, deletes_amount, output_prefix):
    # Creating plot of mutation types and save it
    names = ["inserts", "replaces", "deletes"]
    numbers = [inserts_amount, replaces_amount, deletes_amount]
    plt.bar(names, numbers, align='center', alpha=0.5)
    plt.title('Mutation types')
    plt.ylabel('Appearance')
    plt.xlabel('Types')
    for i in range(len(names)):
        plt.text(x=i, y=numbers[i] + 0.1, s=numbers[i])
    plt.savefig(output_prefix + "-total.png")
    plt.close()


def inserts_diagram(mutations_report, output_prefix):
    # Creating plot of the inserts types and save it
    plt.bar(list(mutations_report.inserts.keys()), list(mutations_report.inserts.values()), align='center', alpha=0.5)
    plt.title('Inserts')
    plt.ylabel('Number Of Inserts')
    plt.xlabel('Bases Inserted')
    for i in range(len(list(mutations_report.inserts.keys()))):
        plt.text(x=i, y=list(mutations_report.inserts.values())[i] + 0.1, s=list(mutations_report.inserts.values())[i])
    plt.savefig(output_prefix + "-inserts.png")
    plt.close()


def replaces_diagram(mutations_report, output_prefix):
    # Creating plot of the replaces types and save it
    plt.bar(list(mutations_report.replaces.keys()), list(mutations_report.replaces.values()), align='center', alpha=0.5)
    plt.title('Replaces')
    plt.ylabel('Number Of Replaces')
    plt.xlabel('Replacement Bases')
    for i in range(len(list(mutations_report.replaces.keys()))):
        plt.text(x=i, y=list(mutations_report.replaces.values())[i] + 0.1,
                 s=list(mutations_report.replaces.values())[i])
    plt.savefig(output_prefix + "-replaces.png")
    plt.close()


def deletes_diagram(mutations_report, output_prefix):
    # Creating plot of the deletes types and save it
    plt.bar(list(mutations_report.deletes.keys()), list(mutations_report.deletes.values()), align='center', alpha=0.5)
    plt.title('Deletes')
    plt.ylabel('Number Of Deletes')
    plt.xlabel('Bases Deleted')
    for i in range(len(list(mutations_report.deletes.keys()))):
        plt.text(x=i, y=list(mutations_report.deletes.values())[i] + 0.1, s=list(mutations_report.deletes.values())[i])
    plt.savefig(output_prefix + "-deletes.png")
    plt.close()


def distance_histogram(array_of_lengths, prefix, title="Histogram"):
    fig = plt.figure()
    plt.hist(array_of_lengths, facecolor='g')
    plt.title(title)
    plt.xlabel('Edit Distance')
    plt.ylabel('Quantity')
    # plt.grid(True)
    fig.savefig(prefix+'.png')



