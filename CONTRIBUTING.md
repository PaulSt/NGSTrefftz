# Contributing to `NGSTrefftz`

🎉**Thanks for taking the time to contribute!** 🎉

To get an overview of the project, check out the [README](README.md).

The [issue tracker](https://github.com/paulst/ngstrefftz/issues)
is the preferred channel for [bug reports](#bugs), [features
requests](#features) and [submitting pull requests](#pull-requests). 

For personal support requests you can reach me at [p.stocker@math.uni-goettingen.de](mailto:p.stocker@math.uni-goettingen.de).

## Bug reports

A bug is a _demonstrable problem_ that is caused by the code in the repository.
Bug reports are extremely helpful - thank you!

Guidelines for bug reports:

1. **Check if the issue has been fixed**: try to reproduce it using the latest `main` or development branch in the repository.

2. **Use the GitHub issue search**: check if the issue has already been reported.

3. **Isolate the problem**: Create a minimal example showing the problem.

4. **Open an issue**: Using the [issue tracker](https://github.com/paulst/ngstrefftz/issues), describe the expected outcome and report the OS, the compiler, and NGSolve version you are using.

## Pull requests

Pull requests - patches, improvements, new features - are a fantastic
help. They should remain focused in scope and avoid containing unrelated commits.
**Please ask first** before embarking on any significant pull request.

Tips on opening a pull request:

1. [Fork](http://help.github.com/fork-a-repo/) the project.

2. Create a branch and implement your feature.
   ```bash
   git checkout -b <your-feature-name>
   ```

3. Run the test suite by calling 
   ```bash
   make test
   ```
   in your build directory. Consider adding new tests for you feature - have a look in the test folder.

4. Once the implementation is done, use Git's
   [interactive rebase](https://help.github.com/articles/interactive-rebase)
   feature to tidy up your commits.
   ```bash
   git rebase --interactive --fork-point main <your-feature-name> 
   ```

5. Push your topic branch up to your fork and [open a Pull Request](https://help.github.com/articles/using-pull-requests/).

**IMPORTANT**: By submitting a patch, you agree to allow the project owners to license your work under the terms of the *LGPL License*.

## Code format

A code style file is provided in the repo to be used with `clang-format`. All c++ files should be formated using this style, to quickly apply it to all files use
   ```bash
   find src/ -iname *.hpp -o -iname *.cpp | xargs clang-format -style='file:.clang-format' -i
   ```
