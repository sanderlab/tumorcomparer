

### GSoC 2021 Reports
<table>
  {% for post in site.reports reversed %}
    {% if post.tags contains "gsoc" %}
      <tr>
        <td><small style="width:5em">{{post.date | date: "%b %d, %y"}}</small></td>
        <td><a href="{{ post.url | remove_first:'/' }}">{{ post.title }}</a></td>
      </tr>
    {% endif %}
  {% endfor %}
</table>


### Contributers
<ul>
  {% for member in site.data.contributers %}
      <li>
        <a title="{{ member.bio }}" target="_blank" href="https://github.com/{{member.github}}">{{ member.name }} ({{ member.position }})</a>
      </li>
  {% endfor %}
</ul>
