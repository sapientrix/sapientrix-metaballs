# Sapientrix Metaballs — Logo generator

## Local Dev

```
npm i
npm run dev
```

## Deploy to GitHub Pages

This project uses Vite. A ready-to-go GitHub Actions workflow is included to build and deploy to GitHub Pages.

Steps:

1) Push to `main` on GitHub
- Ensure your repo is on GitHub and your default branch is `main`.
- The workflow at `.github/workflows/deploy.yml` runs on every push to `main`.

2) Enable Pages in GitHub
- Go to `Settings → Pages` in your repository.
- Under "Build and deployment", set Source to "GitHub Actions".

3) Wait for the workflow to finish
- It builds with the correct Vite `base` for your repo name and deploys the `dist` folder to Pages.
- The URL will be available in the Actions run summary and under `Settings → Pages`.

Notes:
- If your repository is a user/organization site (e.g. `yourname.github.io`), the workflow automatically uses `/` as the base.
- For a project site (e.g. `yourname.github.io/sapientrix-metaballs`), it automatically sets the base to `/${repo}/`.

Manual trigger:
- You can also run the "Deploy to GitHub Pages" workflow manually via the Actions tab (`Run workflow`).
