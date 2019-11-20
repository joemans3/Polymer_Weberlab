# Polymer


## Getting Started
TODO
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites
TODO
What things you need to install the software and how to install them

```
Boost library, standard library, c++11 compiler, libconfig
```

### Installing
TODO
Compile with command
```
g++ -std=c++11 -lconfig++ temp_pf.cpp -o your_file_name
```
Run your executible with alone and have the program prompt you to input parameters or with an init file in the format given. For latter use command (after compiling):
```
./your_file_name init_file.cfg

```
Comes with compiled OSX v.10.14.4. Have not tested with windows but linux works.


Alternative Route:
You can also let Cmake/make compile for you.

After cloning this repository create a /Polymer_Weberlab/build/ directory in the main directory and change directory to it in Terminal. 

Type in the following to create make files:

```
$ cmake ..
```

Followed by:

```
$ make
$ make clean
```

This will create an executable in the build directory called Polymer. You can interact with this exe using the above commands or use a custom script which creates logs and also iterates over multiple simulations. This is the file ```script_polymer``` (bash) in the /Polymer_Weberlab/tools/ directory. The output of this file will be in the same directory in terms of both log files and data files from Polymer. 

The script is currently set up to iterate over multiple temperatures for a set simulation structure but the user can easily change this to anything they like using the current code as template (and a little bash knowhow). TODO: create other custom bash scripts.

<!---
A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds
--->
## Contributing
TODO
<!---
Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.
--->
## Versioning
TODO
<!---
We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 
--->
## Authors

* **Baljyot Parmar**

See also the list of [contributors](https://github.com/joemans3/C_Polymer/contributors) who participated in this project.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details

## Acknowledgments
[Stack Overflow](https://stackoverflow.com) and other mentioned users in specific code blocks.  
