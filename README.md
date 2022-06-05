# plotterdays
plot plot.

# Tools in use
- `vsketch` -- Python library + CLI to create SVG sketches for plotters
- `vpype` -- CLI to help morp SVGs for plotting
- `axicli` -- CLI to interact with the Axidraw plotter

# Development

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


### Common Issues

Sketch redrawing wasn't working for me with Vim / Ubuntu. After tracing throught the vsketch_cli codebase I noticed that it used a `pip` package called `watchfiles` (https://watchfiles.helpmanual.io/#how-watchfiles-works) which in turn uses a Rust package called `notify-rs` (https://github.com/notify-rs/notify). After watching a particular file and inspecting the events coming back I noticed that the system was seing "modified" events but those events weren't being propagated back to the
notify-rs package and in turn the `watchfiles` package. This meant the Qt event handler defined in the `vsketch_cli` wasn't able to trigger redraws.

After searching for a while I found this bug that addresses the issue: https://github.com/notify-rs/notify/issues/394. Ultimately the problem lies with Vim's `backupcopy` option and what it's set to. By default Vim creates a new file and writes to it (this results in a different inode being created and thus tricks the notifiers). See this comment: https://github.com/notify-rs/notify/issues/247#issuecomment-643818591 
