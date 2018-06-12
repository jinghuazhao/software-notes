# Jannovar

From the GitHub repository, it is seen to use `project object model` (POM), an XML representation of a Maven project held in a file named `pom.xml`. We therefore install `maven` first,
```bash
sudo apt install maven
```
The installation then proceeds as follows,
```bash
git clone https://github.com/charite/jannovar
cd jannovar
mvn compile
mvn test
```
