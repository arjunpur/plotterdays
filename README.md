# plotterdays
plot plot.

# Tools in use
- `vsketch` -- Python library + CLI to create SVG sketches for plotters
- `vpype` -- CLI to help morp SVGs for plotting
- `axicli` -- CLI to interact with the Axidraw plotter

# Development

*Raise the pen up / down*

```
axicli -m toggle
```

*Creating a new plot*

```
vsk init -l -p=letter hello_world
```

*Running a sketch*

```
vsk run sketch hello_world

cd plots/hello_world && vsk run

vsk run sketch_hello_world.py
```

*Saving a sketch*

```
vsk save
```

*Using axicli*

Sending a sketch to Axidraw
```
axicli output/hello_world.svg
```

Enable / Disable the motors:
```
axicli -m manual --manual_cmd=[enable_xy|disable_xy
```

*Drawing with layers*

In code you can use `vsk.stroke(<layer_number>)` to draw at a particular layer. You can then use
the following commands to draw the layers you want:

```
axicli -m=layers --layer=1 --speed_pendown=15 spirals.svg
axicli -m=layers --layer=2 --speed_pendown=15 spirals.svg
```

*Keeping the paper aligned with the plotter*
I've been using a ruler to measure the distance between the plotter and the paper / MDF board at many points along the width of the paper (parallel to the plotter). The distance between the plotter and the paper should be the same if they're aligned.


### Common Issues / FAQ / Knowledge

- Sketch redrawing wasn't working for me with Vim / Ubuntu. After tracing throught the vsketch_cli codebase I noticed that it used a `pip` package called `watchfiles` (https://watchfiles.helpmanual.io/#how-watchfiles-works) which in turn uses a Rust package called `notify-rs` (https://github.com/notify-rs/notify). After watching a particular file and inspecting the events coming back I noticed that the system was seing "modified" events but those events weren't being propagated back to the
notify-rs package and in turn the `watchfiles` package. This meant the Qt event handler defined in the `vsketch_cli` wasn't able to trigger redraws.

- After searching for a while I found this bug that addresses the issue: https://github.com/notify-rs/notify/issues/394. Ultimately the problem lies with Vim's `backupcopy` option and what it's set to. By default Vim creates a new file and writes to it (this results in a different inode being created and thus tricks the notifiers). See this comment: https://github.com/notify-rs/notify/issues/247#issuecomment-643818591 

- In order for your python project to be able to resolve imports of common libraries, you need to make sure the project directory is in the PYTHONPATH. I have the following:
```
export PYTHONPATH="$HOME/Projects/:$PYTHONPATH"
```

- `vsketch` is installed as a standalone application using `pipx`. There are some complexities with virtualenv because `vsketch` installs packages in it's own virtualenv directory, whereas packages that I install are in a different virtualenv directory

- `axicli` gets installed in ~/.pyenv/bin/

- Make sure you set `MYPYPATH` to your project directory so that it can find your modules

#### Upgrading Python Version
- https://stackoverflow.com/questions/44692668/python-how-can-i-update-python-version-in-pyenv-virtual-environment
- When doing this, you'll need to manually install the axidrawcli: https://axidraw.com/doc/cli_api/#introduction. Some versions will collide so remove the axidraw dependencies from the requirements.txt, install everything, and then install axidraw.


### Resources

- AxiDraw user guide: https://cdn.evilmadscientist.com/dl/ad/public/AxiDraw_Guide_v550.pdf


### dependencies

For rendering LaTex to SVGs:
- Latex (Tex Live for Ubuntu) (`sudo apt update && sudo apt install texlive`)
- LuaTex

- vsketch
- vpype
- Inkscape (?)
- axicli

- Python3.10+
