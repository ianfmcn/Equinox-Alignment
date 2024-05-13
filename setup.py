setup(
    name="equinox",
    version=0.1,
    description="CSE 185 Project",
    author="Natsuki Romero, Leah Kim, Ian McNellis",
    author_email=""
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "equinox=equinox.equinox:main"
        ],
    },
)
