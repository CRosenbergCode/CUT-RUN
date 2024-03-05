## Templates for Alpine


These scripts contain useful SBATCH settings to defaults. Also included is setting TMP environmental variable
to scratch space on the runtime node (versus alpine scratch), which can be used as temporary disk space by commands in your script.


### Array jobs

Sometimes the easiest way to parallelize a workflow is to run the same script multiple times with different arguments, or on different data.
The `sbatch --array` provides a shortcut for doing so, as long as your script makes use of the `$SLURM_ARRAY_TASK_ID` environmental variable
to distinguish the arguments you wish to vary between jobs.

#### Parallelize over the command line 

For example, running the following:

```
sbatch --array=0-4 array-over-cmdline-args.sbatch apples oranges grapes bananas cherries
```

will submit a separate job 5 times, but the script uses the `$SLURM_ARRAY_TASK_ID` to choose one of the command line arguments (numbered 0 through 4), and
pass it to a command.

| submitted command line | `$SLURM_ARRAY_TASK_ID` | value of `arg` in script |
-------------------------|-----------------------|---------------------------
| array-over-cmdline-args.sbatch apples oranges grapes bananas cherries | 0 | apples  |
| array-over-cmdline-args.sbatch apples oranges grapes bananas cherries | 1 | oranges |
| array-over-cmdline-args.sbatch apples oranges grapes bananas cherries | 2 | grapes  |
| array-over-cmdline-args.sbatch apples oranges grapes bananas cherries | 3 | bananas |
| array-over-cmdline-args.sbatch apples oranges grapes bananas cherries | 4 | cherries|

Look at the script `array-over-cmdline-args.sbatch` to see how `arg` is extracted from the command line arguments using `$SLURM_ARRAY_TASK_ID`.

#### Parallelize over lines in an input file

Sometimes it is convenient to specify sets of arguments for each run in a file.  Using a `while` loop, we can iterate over lines in a file and check if
our iteration number matches `$SLURM_ARRAY_TASK_ID`. We execute the command if it does.

The following parallelizes execution of a command once per line of input in `example-metadata.txt`:

`sbatch --array=0-3 array-loop-over-file-template.sbatch`

| submitted command line | `$SLURM_ARRAY_TASK_ID` | value of `fields` in script |
-------------------------|-----------------------|---------------------------
| array-over-cmdline-args.sbatch | 0 |  andersen    darling |
| array-over-cmdline-args.sbatch | 1 |  david       duchovny |
| array-over-cmdline-args.sbatch | 2 |  metropolis  hastings |
| array-over-cmdline-args.sbatch | 3 |  forward     backward |

See the source code of `array-loop-over-file-template.sbatch` for specifics on looping over a file and setting its contents to an array.

Note: if you have a very large file, it might be more performant to use a utility like `awk` or `sed` to extract the line of interest, rather than a loop.
