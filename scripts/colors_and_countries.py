country_list = [
    "France",
    "United Kingdom",
    "Austria",
    "Netherlands",
    "Norway",
    "Spain",
    "Belgium",
    "Switzerland",
    "Germany",
    "USA",
    "Ireland",
    "Italy",
    "New Zealand",
    "Sweden",
    "Australia",
    "Singapore",
    "Denmark",
    "Hong Kong",
    "Portugal",
    "Luxembourg",
    "Lithuania",
    "Iceland",
    "Czech Republic",
    "Israel",
    "Finland",
    "Poland",
]

colors = [
    "#a6cee3",
    "#1f78b4",
    "#ffff00",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#b2df8a",
    "#b15928",
    "#000000",
    "#9900cc",
    "#ffcc00",
    "#a6a6a6",
    "#0099ff",
    "#d9d9d9",
    "#990000",
    "#c2c2a3",
    "#ff66ff",
    "#0000ff",
    "#006600",
    "#cc00cc",
    "#993300",
    "#336699",
]

# separate out Scotland, England, NI, Wales...
uk_countries = ["Scotland", "England", "Wales", "Northern Ireland"]
countries_and_uk_list = [
    x for x in country_list if x != "United Kingdom"
] + uk_countries
all_countries = country_list + uk_countries

# For the main cluster, ran out of colors... try this hack.
# palette = sns.color_palette("tab10", round(len(country_list)/2))
linestyles = ["-", "-.", "--", ":"]

country_styles = {
    country: {"c": colors[i], "ls": linestyles[i // len(colors)]}
    for i, country in enumerate(country_list)
}

country_styles.update(
    {
        "Scotland": {"c": country_styles["United Kingdom"]["c"], "ls": "--"},
        "England": {"c": country_styles["United Kingdom"]["c"], "ls": "-"},
        "Wales": {"c": country_styles["United Kingdom"]["c"], "ls": "-."},
        "Northern Ireland": {"c": country_styles["United Kingdom"]["c"], "ls": ":"},
    }
)

popsizes = {
    "Spain": 46723749,
    "United Kingdom": 66488991,
    "Switzerland": 8516543,
    "Netherlands": 17231017,
    "Norway": 5314336,
    "Germany": 82927922,
    "France": 66987244,
    "Ireland": 4853506,
    "Belgium": 11422068,
    "Denmark": 5797446,
    "Scotland": 5469000,
    "Wales": 3136000,
}

colors_2 = [
    "#a6cee3",
    "#1f78b4",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#b2df8a",
    "#b15928",
    "#9900cc",
    "#990000",
    "#ff66ff",
    "#0000ff",
    "#a6a6a6",
    "#000000",
    "#004d00",
    "#ffcc00",
    "#fdbf6f",
    "#0099ff",
]

country_list_1 = [
    "Czech Republic",  # 1   # "#a6cee3",
    "United Kingdom",  # 2   # "#1f78b4",
    "Netherlands",  # 3      # "#33a02c",
    "Norway",  # 4           # "#fb9a99",
    "India",  # 5            # "#e31a1c",
    "Switzerland",  # 6      # "#ff7f00",
    "Germany",  # 7          # "#cab2d6",
    "USA",  # 8              # "#6a3d9a",
    "Ireland",  # 9          # "#b2df8a",
    "Italy",  # 10           # "#b15928",
    "Russia",  # 11          # "#9900cc",
    "Portugal",  # 12        # "#990000",
    "Lithuania",  # 13       # "#ff66ff",
    "Spain",  # 14           # "#0000ff",
    "Singapore",  # 15       # "#a6a6a6",
    "Thailand",  # 16        # "#000000",
    "France",  # 17          # "#004d00",
    "Australia",  # 18       # "#ffcc00",
    "Belgium",  # 19         # "#fdbf6f",
    "Denmark",  # 20         # "#0099ff",
]

country_list_2 = [
    "Brazil",  # 1
    "Slovakia",  # 2
    "Botswana",  # 3
    "Turkey",  # 4
    "South Africa",  # 5
    "Curacao",  # 6
    "Sweden", # 7
    "Poland",  # 8
    "Romania",  # 9
    "Canada",  # 10
    "Iceland",  # 11
    "Colombia",  # 12
    "Ghana",  # 13
    "Slovenia",  # 14
    "Latvia",  # 15
    "Nigeria",  # 16
    "Aruba",  # 17
    "Austria",  # 18
    "Luxembourg",  # 19
    "Israel",  # 20
]

country_list_3 = [
    "Finland",  # 1
    "Japan",  # 2
    "Croatia",  # 3
    "Bangladesh",  # 4
    "Mexico",  # 5
    "South Korea",  # 6
    "Malaysia",  # 7
    "Bulgaria",  # 8
    "Indonesia",  # 9
    "Estonia",  # 10
    "Greece",  # 11
    "Chile",  # 12
    "Qatar",    # 13
    "Sri Lanka",   # 14
    "Egypt",    # 15
    "Philippines",     # 16
    "Argentina", # 17
    "Cambodia",       # 18  
    "Mozambique",    #19
    "Zimbabwe",    #20
]

country_list_4 = [
    "Kenya",                #1
    "Ecuador",              #2
    "Malawi",               #3
    "Sint Maarten",         #4
    "Angola",               #5
    "Trinidad and Tobago",  #6
    "Peru",                 #7
    "North Macedonia",      #8
]

country_styles_1 = {
    country: {"c": colors_2[i], "ls": "-"}  # linestyles[i//len(colors)]}
    for i, country in enumerate(country_list_1)
}

country_styles_2 = {
    country: {"c": colors_2[i], "ls": "--"}  # linestyles[i//len(colors)]}
    for i, country in enumerate(country_list_2)
}

country_styles_3 = {
    country: {"c": colors_2[i], "ls": "-."}  # linestyles[i//len(colors)]}
    for i, country in enumerate(country_list_3)
}

country_styles_4 = {
    country: {"c": colors_2[i], "ls": "dotted"}  # linestyles[i//len(colors)]}
    for i, country in enumerate(country_list_4)
}

country_styles_all = {**country_styles_1, **country_styles_2, **country_styles_3, **country_styles_4}
