# Uploading this project to GitHub (GitHub Desktop)

The supplied folder is **not yet published on GitHub**. For a first upload:

1. Extract the downloaded ZIP. Keep the top-level `epilepsy-phd-analysis-code` folder intact.
2. Install and sign in to **GitHub Desktop** (https://desktop.github.com/).
3. Choose **File → Add local repository** and select the extracted top-level folder. If prompted that it is not yet a Git repository, select the option to **create a repository here**. If your Desktop version does not offer that option, choose **File → New repository**, and select the existing project folder/location without creating a duplicate nested folder.
4. Check the changed-files list. You should see code, small CSV/TSV files and Markdown documentation, **not** raw GWAS files, `.rds` files or original HPC `.log` files.
5. Enter a commit summary such as **Initial verified Chapter 3 standalone LDSC code**, then click **Commit to main**.
6. Click **Publish repository**. Use the name `epilepsy-phd-analysis-code`; choose **Private** until you have reviewed it for release.
7. Open your new repository in a browser and confirm that the README and Chapter 3 directories are visible.

Do **not** upload the ZIP as one file within a GitHub repository. Publish the **contents of the extracted folder** using GitHub Desktop.

For future chapters, add their verified files into the same local repository, commit the changes, and push from GitHub Desktop. When your thesis is ready to hyperlink the repository, switch visibility to public only after an institutional/data-sharing check and add a citation and software license if appropriate.

If you already have a preferred repository name, rename the repository when publishing; the local folder name need not match exactly.
