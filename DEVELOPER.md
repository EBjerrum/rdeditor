# Developer Guide

## Installation

### Prerequisites

rdkit, numpy, pyside2 are installed by default.
ruff is installed using the dev tag as described below

### Installation Steps

1. Clone the repository:

   ```bash
   git clone https://github.com/EBjerrum/rdeditor.git
   ```

2. Navigate to the project directory:

   ```bash
   cd rdeditor
   ```

3. Install rdeditor in editable mode to enable developer modifications:

   ```bash
   pip install -e .[dev]
   ```

4. Optionally, set up your preferred code editor (e.g., VS Code) to format code on save using ruff.

## Automatic Code Formatting with Ruff

The `rdeditor` project utilizes `ruff` as both a linter and code formatter to ensure consistent code quality and formatting standards. You can automate the code formatting process by setting up a pre-commit hook in your local Git repository or configuring VS Code to auto-format code on save using the provided `ruff.toml` specifications.

### Creating a Pre-Commit Hook

You can set up a pre-commit hook in your local Git repository to automatically format code using `ruff` before each commit. Here's how to do it:

1. Navigate to the `.git/hooks` directory in your repository.

2. Create a new file named `pre-commit` if it doesn't already exist.

3. Open the `pre-commit` file in a text editor and add the following content:

   ```bash
   #!/bin/bash

   # Run the ruff formatter on staged changes
   ruff format $(git diff --cached --name-only | grep '\.py$')

   # Stage the formatted changes
   git add $(git diff --cached --name-only)
   ```

4.

5. Save the file and make it executable by running the following command in your terminal:

   ```bash
   chmod +x .git/hooks/pre-commit
   ```

Now, each time you attempt to commit changes, the `pre-commit` hook will run the `ruff` formatter on your staged changes, ensuring consistent formatting before the commit is finalized.

### Configuring VS Code for Auto-Formatting on Save

If you prefer to use VS Code, you can configure it to automatically format code using the provided `ruff.toml` specifications on save. Here's how to do it:

1. Open VS Code and navigate to the settings by clicking on the gear icon in the bottom left corner or by pressing `Ctrl + ,`.

2. In the search bar at the top, type "format on save" to find the setting.

3. Check the box next to "Editor: Format On Save" to enable auto-formatting on save.

4. Next, click on "Extensions" in the sidebar and search for "ruff" in the search bar.

5. Install the ruff extension.

6. Press `ctrl-shift-p` or select `format with...` from the right click menu in a python file. Select the option `configure default formatter...` and choose ruff.

With these settings configured, VS Code will automatically format your Python code according to the specifications provided in the `ruff.toml` file each time you save a file.

## Usage

[Include usage instructions here if different from README.md]

## Contributing

## Checking code with ruff

Code submitted to github master branch are checked with ruff via GitHub Actions. It is thus advisable to check yourself before a pull request is made. This can be done with:

`ruff check` which will inspect the code and print a list of issues.

Sometimes they can be safely fixed with

`ruff check --fix` and even `ruff check --fix --unsafe-fixes` but otherwise they need to be inspected and mitigated.

### Development Environment Setup

1. Follow the installation steps in the INSTALL section.
2. Set up your preferred code editor (e.g., VS Code) to format code on save.
3. Optionally, configure a pre-commit hook locally to ensure code consistency before committing changes.

### Branching Strategy

- Use meaningful branch names for new features, bug fixes, or enhancements (e.g., feature/add-new-tool, fix/issue-123).
- Create feature branches from the main branch and submit pull requests for review.
- Review and address any feedback from maintainers before merging changes.

### Code Guidelines

- Follow the project's coding standards and style guide.
- Write clear, concise, and well-documented code.
- Include relevant comments and docstrings to explain complex logic or functionality.
- Write meaningful commit messages that describe the purpose of the changes.

## Testing

[Include testing instructions and guidelines here if applicable]

## Troubleshooting

[Include troubleshooting tips and common issues here]

## Additional Resources

[Include links to relevant documentation, tutorials, or external resources]

## Contact

[Provide contact information for project maintainers or contributors]

## Changelog

[Include a summary of recent changes, improvements, and bug fixes]
