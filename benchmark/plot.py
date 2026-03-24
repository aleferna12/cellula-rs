import json
import plotly.express as px

with open("benchmark/results.json") as f:
    results = json.load(f)["results"]
stats = {
    "max_pop": [],
    "size": [],
    "mean": [],
    "stddev": [],
    "framework": []
}
for run in results:
    for param in run["parameters"]:
        stats[param].append(run["parameters"][param])
    stats["mean"].append(run["mean"])
    stats["stddev"].append(run["stddev"])
    command = run["command"]
    if "target/fastest" in command:
        command = "cellulars"
    elif "morpheus" in command:
        command = "morpheus"
    stats["framework"].append(command)

fig = px.line(
    x=stats["max_pop"],
    y=stats["mean"],
    error_y=stats["stddev"],
    line_dash=stats["size"],
    color=stats["framework"],
    width=600,
    height=400,
    labels={"color": "framework", "line_dash": "lattice size"},
).update_layout(
    template="plotly_white",
    xaxis_type="log",
    xaxis_title="nr. of cells",
    yaxis_title="execution time",
)
fig.write_html("benchmark/results.html")
fig.write_image("benchmark/results.svg")