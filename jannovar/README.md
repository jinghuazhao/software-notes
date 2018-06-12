# Jannovar

From the GitHub repository, it is seen to use `project object model` (POM), an XML representation of a Maven project held in a file named `pom.xml`. We therefore install `maven` first,
```bash
sudo apt install maven
```
The installation then proceeds as follows,
```bash
git clone https://github.com/charite/jannovar
cd jannovar
mvn package
```
Other tasks such as compile, test, etc. are also possible.

It is handy to use symbolic link, i.e.,
```bash
ln -s /home/jhz22/D/genetics/jannovar/jannovar-cli/target/jannovar-cli-0.24.jar $HOME/bin/Jannovar.jar
java -jar $HOME/bin/Jannovar.jar db-list
java -jar Jannovar.jar download -d hg19/refseq
```
