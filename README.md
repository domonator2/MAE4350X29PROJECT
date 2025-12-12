New Age Aerospace - X-29 Reverse Engineering Project (MAE 4350)
Welcome to the official code repository for the New Age Aerospace team. This project contains all MATLAB code and supporting files for the reverse engineering of the Grumman X-29.

The primary goal of this repository is to ensure organized, conflict-free collaboration. All team members are required to read and follow the workflow outlined below.

📂 Repository Structure
All code is organized into folders corresponding to the project disciplines. The Chief Engineer (Synthesis) is the only person who will merge code into the main_simulation.m file.

/Aerodynamics/: Code and data related to aerodynamic analysis. (Lead: Andrea)

/Cost/: Code for cost modeling and analysis. (Lead: Tan)

/Geometry/: Scripts for generating and analyzing vehicle geometry (OpenVSP files, etc.). (Lead: Jayden)

/Propulsion/: Engine models and propulsion system analysis. (Lead: Carlos)

/Stability/: Stability and control analysis scripts. (Lead: Ethan)

/Structures_WB/: Scripts for weights, balance, and structural analysis. (Lead: Samantha)

/Trajectory_Performance/: Mission profile and performance simulation code. (Lead: Damian)

/Synthesis/: Integration scripts and top-level functions managed by the Chief Engineer. (Lead: Dominic)

/Common/: Any shared functions, constants, or data files that are used by more than one discipline.

main_simulation.m: (DO NOT EDIT DIRECTLY) This is the main, top-level script that runs the entire simulation. It will be managed only by the Synthesis lead.

⚙️ The Workflow: How We Work Together
To prevent code conflicts and lost work, we will use a standard "branching" workflow. No one should ever work directly on the main branch.

Your 5-Step Process for Adding Code:
Create a New Branch: Before you start writing code, create your own branch.

In MATLAB or on the GitHub website, create a new branch from main.

Name it descriptively: yourname-feature. For example: andrea-drag-calculation.

Write Your Code: Work on your MATLAB scripts within your branch. You can make as many changes as you need.

Commit Your Changes: Save your work by "committing" it.

Each commit should be a small, logical change.

Use a clear commit message, e.g., "Added function to calculate wing lift coefficient."

Open a Pull Request (PR): When your feature is complete and tested, open a "Pull Request".

This is a formal request to merge your branch into the main branch.

On GitHub, go to the "Pull Requests" tab and click "New pull request".

Set the "base" branch to main and the "compare" branch to your feature branch.

Assign the Chief Engineer (Dominic) as the reviewer.

==TEST==

Review and Merge: The Chief Engineer will review your code to ensure it integrates properly with the main simulation. Once approved, your code will be merged.

IMPORTANT: Before starting new work, always update your local main branch to get the latest changes from the team!
