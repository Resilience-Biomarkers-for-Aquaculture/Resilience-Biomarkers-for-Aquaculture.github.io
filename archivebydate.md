---
layout: page
title: Post by Month
permalink: /monthview/
active: archivebydate
sitemap: false
---

<div id="index">
{% for post in site.posts %}
  {% capture year %}{{ post.date | date: '%Y' }}{% endcapture %}
  {% capture nyear %}{{ post.next.date | date: '%Y' }}{% endcapture %}
  {% unless post.next %}
    <h2 class="archivetitletopbottom">{{ post.date | date: '%Y' }}</h2>
  {% else %}
    {% if year != nyear %}
      {% if forloop.index != 1 %}</ul><!-- year -->{% endif %}
      <h2 class="archivetitletopbottom">{{ post.date | date: '%Y' }}</h2>
    {% endif %}
  {% endunless %}

  {% capture month %}{{ post.date | date: '%m%Y' }}{% endcapture %}
  {% capture nmonth %}{{ post.next.date | date: '%m%Y' }}{% endcapture %}
  {% if month != nmonth %}
    {% if forloop.index != 1 && year == nyear %}</ul><!-- month -->{% endif %}
    <h2 class="archivetitle">{{ post.date | date: '%B %Y' }}</h2><ul>
  {% endif %}

  {% if post.link %}
    <h3 class="link-post">
      <a href="{{ site.baseurl }}{{ post.url }}" title="{{ post.title }}">{{ post.title }}</a>
      <a href="{{ post.link }}" target="_blank" title="{{ post.title }}">
        <i class="fa fa-link"></i>
      </a>
    </h3>
  {% else %}
    <li>
      <a href="{{ post.url | prepend: site.baseurl }}">{{ post.title }}</a>
      {% if post.author %} • {{ post.author }}{% endif %}
      {% if post.date %} • {{ post.date | date: "%B %e, %Y" }}{% endif %}
    </li>
  {% endif %}
{% endfor %}
</ul>

</div>