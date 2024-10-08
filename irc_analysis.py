import re
import argparse
import os
import matplotlib.pyplot as plt

def parse_irc_data(qchem_output):
    """
    Parses the IRC data from the Q-Chem output.

    Args:
        qchem_output (str): The content of the Q-Chem output file.

    Returns:
        tuple: A tuple containing two lists:
            - irc_steps: A list of IRC step data. Each element is a dictionary containing:
                - 'step': IRC step number
                - 'energy': Total energy at the step
                - 'coordinates': A list of tuples representing atom coordinates in XYZ format
            - energies: A list of energies for the plot
    """
    irc_steps = []
    energies = []

    lines = qchem_output.splitlines()
    step_data = None

    for index, line in enumerate(lines):
        if "IRC step #" in line:
            num_atoms = int(lines[index-1].split()[0])
            processed_line = line.split()
            irc_step_number = processed_line[3].strip(',')
            energy = processed_line[6]
            if step_data:
                irc_steps.append(step_data)
            step_data = {
                'step': int(irc_step_number),
                'energy': float(energy),
                'coordinates': []
            }
            energies.append(float(energy))

            for i in range(num_atoms):
                step_data['coordinates'].append((
                    str(lines[index+1+i].split()[0]),
                    float(lines[index+1+i].split()[1]),
                    float(lines[index+1+i].split()[2]),
                    float(lines[index+1+i].split()[3])
                ))

    if step_data:
        irc_steps.append(step_data)

    # Separate forward and backward reactions
    forward_steps = [step for step in irc_steps if step['step'] > 0]
    backward_steps = [step for step in irc_steps if step['step'] < 0]

    # Reverse the forward reaction steps to start at the final point
    forward_steps = forward_steps[::-1]

    # Concatenate both reactions
    irc_steps = forward_steps + backward_steps

    # Combine energies
    energies = [step['energy'] for step in irc_steps]

    return irc_steps, energies


def write_xyz_trajectory(filename, irc_steps):
    """
    Writes an XYZ trajectory file based on the parsed IRC steps.

    Args:
        filename (str): The output XYZ file name.
        irc_steps (list): The parsed IRC steps data.
    """
    print("Writing xyz trajectories")
    with open(filename, 'w') as f:
        for step in irc_steps:
            f.write(f"{len(step['coordinates'])}\n")
            f.write(f"IRC step {step['step']}, Energy {step['energy']}\n")
            for atom in step['coordinates']:
                f.write(f"{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n")


def get_relative_energies(energies):
    """
    Find minimum energy, then uses this energy as the zero point, then converts
    all energies to kcal/mol

    Args:
        energies (list): The parsed energy data.
    """
    relative_energies = [None] * len(energies)
    for i in range(len(energies)):
        relative_energies[i] = (energies[i] - min(energies))*627.509

    return relative_energies


def plot_irc_trajectory(energies, output_plot):
    """
    Plots the IRC energies.

    Args:
        energies (list): The list of energies.
        output_plot (str): The output file name for the plot.
    """
    print("Plotting...")
    plt.plot(range(1, len(energies) + 1), energies, marker='o')
    plt.xlabel('IRC Step')
    plt.ylabel('∆E (kcal/mol)')
    plt.title('IRC Energy Profile')
    plt.grid(True)
    plt.savefig(output_plot)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Parse a Q-Chem output file and generate an XYZ trajectory and energy plot.")
    parser.add_argument('input_file', type=str, help="Input Q-Chem output file")

    args = parser.parse_args()
    input_basename = os.path.splitext(os.path.basename(args.input_file))[0]
    xyz_output = f"{input_basename}.xyz"
    plot_output = f"{input_basename}.png"

    # String that signals completion of IRC calculation
    search_string = "IRC -- convergence criterion reached."

    # Load Q-Chem output file
    with open(args.input_file, 'r') as file:
        qchem_output = file.read()

    # Find the index of the last occurrence of the search string
    last_index = qchem_output.rfind(search_string)

    if last_index == -1:
        print(f'Warning: "{search_string}" not found in the file.')
        # Handle case where search_string is not found
        qchem_output_after_last = ""
    else:
        print("Successful termination of IRC calculation!")
        # Extract everything after the last occurrence of the search string
        qchem_output_after_last = qchem_output[last_index:]

    # Parse the IRC data
    irc_steps, energies = parse_irc_data(qchem_output_after_last)

    # Write the XYZ trajectory
    write_xyz_trajectory(xyz_output, irc_steps)

    # Get the relative energies with respect to the lowest energy, in kcal/mol
    relative_energies = get_relative_energies(energies)

    # Plot the IRC trajectory
    plot_irc_trajectory(relative_energies, plot_output)

    print(f"XYZ trajectory and IRC plot generated successfully as {xyz_output} and {plot_output}.")


if __name__ == '__main__':
    main()
